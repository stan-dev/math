#ifndef STAN_MATH_MIX_FUNCTOR_LAPLACE_LIKELIHOOD_HPP
#define STAN_MATH_MIX_FUNCTOR_LAPLACE_LIKELIHOOD_HPP

// #include <stan/math/mix/laplace/hessian_times_vector.hpp>
#include <stan/math/mix/functor/hessian_block_diag.hpp>
#include <stan/math/prim/fun.hpp>
#include <Eigen/Sparse>

namespace stan {
namespace math {
/**
 * functions to compute the log density, first, second,
 * and third-order derivatives for a likelihoood specified by the user.
 */
namespace laplace_likelihood {
namespace internal {
/**
 * @tparam F
 * @tparam Theta
 * @tparam Eta
 * @param f
 * @param theta
 * @param eta
 * @param args
 */
template <typename F, typename Theta, typename Eta, typename... Args,
          require_eigen_vector_t<Theta>* = nullptr,
          require_eigen_t<Eta>* = nullptr>
inline auto log_likelihood(F&& f, const Theta& theta, const Eta& eta,
                           Args&&... args) {
  return f(theta, eta, args...);
}

/**
 * @tparam F
 * @tparam Theta
 * @tparam Eta
 * @tparam Args
 * @param f
 * @param theta
 * @param eta
 * @param gradient
 * @param hessian_block_size
 * @param args
 */
template <typename F, typename Theta, typename Eta, typename... Args,
          require_eigen_vector_t<Theta>* = nullptr,
          require_eigen_t<Eta>* = nullptr>
inline Eigen::SparseMatrix<double> diff(F&& f, const Theta& theta,
                                        const Eta& eta,
                                        Eigen::VectorXd& gradient,
                                        const Eigen::Index hessian_block_size,
                                        Args&&... args) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  const Eigen::Index theta_size = theta.size();
  const Eigen::Index eta_size = eta.size();
  {
    nested_rev_autodiff nested;
    Matrix<var, Dynamic, 1> theta_var = theta;
    Matrix<var, Eta::RowsAtCompileTime, Eta::ColsAtCompileTime> eta_var = eta;

    var f_var = f(theta_var, eta_var, args...);
    grad(f_var.vi_);
    gradient.resize(theta_size + eta_size);
    for (Eigen::Index i = 0; i < theta_size; i++) {
      gradient(i) = theta_var(i).adj();
    }
    for (Eigen::Index i = 0; i < eta_size; i++) {
      gradient(theta_size + i) = eta_var(i).adj();
    }
  }
  if (hessian_block_size == 1) {
    Eigen::VectorXd v = Eigen::VectorXd::Ones(theta_size);
    Eigen::VectorXd hessian_v = hessian_times_vector(f, theta, eta, v, args...);
    Eigen::SparseMatrix<double> hessian_theta(theta_size, theta_size);
    hessian_theta.reserve(Eigen::VectorXi::Constant(theta_size, 1));
    for (Eigen::Index i = 0; i < theta_size; i++) {
      hessian_theta.insert(i, i) = hessian_v(i);
    }
    return hessian_theta;
  } else {
    return hessian_block_diag(f, theta, eta, hessian_block_size, args...);
  }
}

/**
 * @tparam F
 * @tparam Theta
 * @tparam Eta
 * @tparam Args
 * @param f
 * @param theta
 * @param eta
 * @param args
 */
template <typename F, typename Theta, typename Eta, typename... Args,
          require_eigen_vector_t<Theta>* = nullptr,
          require_eigen_t<Eta>* = nullptr>
inline Eigen::VectorXd third_diff(F&& f, const Theta& theta, const Eta& eta,
                                  Args&&... args) {
  nested_rev_autodiff nested;
  const Eigen::Index theta_size = theta.size();
  Eigen::Matrix<var, Eigen::Dynamic, 1> theta_var = theta;
  Eigen::Matrix<fvar<var>, Eigen::Dynamic, 1> theta_fvar(theta_size);
  for (Eigen::Index i = 0; i < theta_size; ++i) {
    theta_fvar(i) = fvar<var>(theta_var(i), 1.0);
  }
  fvar<var> ftheta_fvar = f(theta_fvar, eta, args...);

  Eigen::Matrix<fvar<fvar<var>>, Eigen::Dynamic, 1> theta_ffvar(theta_size);
  for (Eigen::Index i = 0; i < theta_size; ++i) {
    theta_ffvar(i) = fvar<fvar<var>>(theta_fvar(i), 1.0);
  }
  fvar<fvar<var>> ftheta_ffvar = f(theta_ffvar, eta, args...);
  grad(ftheta_ffvar.d_.d_.vi_);
  return theta_var.adj();
}

/**
 * @tparam F
 * @tparam Theta
 * @tparam Eta
 * @tparam Args
 * @param f
 * @param theta
 * @param eta
 * @param A
 * @param hessian_block_size
 * @param args
 */
template <typename F, typename Theta, typename Eta, typename... Args,
          require_eigen_vector_t<Theta>* = nullptr,
          require_eigen_t<Eta>* = nullptr>
inline Eigen::VectorXd compute_s2(F&& f, const Theta& theta, const Eta& eta,
                                  const Eigen::MatrixXd& A,
                                  int hessian_block_size, Args&&... args) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::MatrixXd;
  using Eigen::VectorXd;

  nested_rev_autodiff nested;
  const Eigen::Index theta_size = theta.size();
  const Eigen::Index eta_size = eta.size();
  const Eigen::Index parm_size = theta_size + eta_size;
  Matrix<var, Dynamic, 1> theta_var = theta;
  Matrix<var, Eta::RowsAtCompileTime, Eta::ColsAtCompileTime> eta_var = eta;
  int n_blocks = theta_size / hessian_block_size;
  fvar<fvar<var>> target_ffvar = 0;
  VectorXd v(theta_size);
  VectorXd w(theta_size);
  for (Eigen::Index i = 0; i < hessian_block_size; ++i) {
    v.setZero();
    for (int j = i; j < theta_size; j += hessian_block_size) {
      v(j) = 1;
    }
    Matrix<fvar<var>, Dynamic, 1> theta_fvar(theta_size);
    for (int j = 0; j < theta_size; ++j) {
      theta_fvar(j) = fvar<var>(theta_var(j), v(j));
    }
    w.setZero();
    for (int j = 0; j < n_blocks; ++j) {
      for (int k = 0; k < hessian_block_size; ++k) {
        w(k + j * hessian_block_size)
            = A(k + j * hessian_block_size, i + j * hessian_block_size);
      }
    }
    Matrix<fvar<var>, Eta::RowsAtCompileTime, Eta::ColsAtCompileTime> eta_fvar
        = eta_var.template cast<fvar<var>>();
    Matrix<fvar<fvar<var>>, Dynamic, 1> theta_ffvar(theta_size);
    for (int j = 0; j < theta_size; ++j) {
      theta_ffvar(j) = fvar<fvar<var>>(theta_fvar(j), w(j));
    }
    Matrix<fvar<fvar<var>>, Eta::RowsAtCompileTime, Eta::ColsAtCompileTime>
        eta_ffvar = eta_fvar.template cast<fvar<fvar<var>>>();

    fvar<var> f_fvar = f(theta_fvar, eta_fvar, args...);
    target_ffvar += f(theta_ffvar, eta_ffvar, args...);
  }
  grad(target_ffvar.d_.d_.vi_);
  VectorXd parm_adj(parm_size);
  for (Eigen::Index i = 0; i < theta_size; ++i) {
    parm_adj(i) = theta_var(i).adj();
  }
  for (Eigen::Index i = 0; i < eta_size; ++i) {
    parm_adj(theta_size + i) = eta_var(i).adj();
  }
  return 0.5 * parm_adj;
}

/**
 * @tparam F
 * @tparam Theta
 * @tparam Eta
 * @tparam Args
 * @param f
 * @param v
 * @param theta
 * @param eta
 * @param args
 */
template <typename F, typename V_t, typename Theta, typename Eta,
          typename... Args, require_eigen_vector_t<Theta>* = nullptr,
          require_eigen_t<Eta>* = nullptr>
inline plain_type_t<Eta> diff_eta_implicit(F&& f, const V_t& v,
                                           const Theta& theta, const Eta& eta,
                                           Args&&... args) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::VectorXd;
  if constexpr (Eta::RowsAtCompileTime == 0 && Eta::ColsAtCompileTime == 0) {
    return plain_type_t<Eta>{};
  }
  nested_rev_autodiff nested;
  const Eigen::Index eta_size = eta.size();
  Matrix<var, Eta::RowsAtCompileTime, Eta::ColsAtCompileTime> eta_var = eta;

  // CHECK -- can we avoid declaring theta as fvar<var>?
  // We currently compute derivatives wrt eta, which is not needed.
  const Eigen::Index theta_size = theta.size();
  Matrix<var, Dynamic, 1> theta_var = theta;
  Matrix<fvar<var>, Dynamic, 1> theta_fvar(theta_size);
  for (Eigen::Index i = 0; i < theta_size; i++) {
    theta_fvar(i) = fvar<var>(theta_var(i), v(i));
  }
  Matrix<fvar<var>, Eta::RowsAtCompileTime, Eta::ColsAtCompileTime> eta_fvar
      = eta_var.template cast<fvar<var>>();

  fvar<var> f_fvar = f(theta_fvar, eta_fvar, args...);
  grad(f_fvar.d_.vi_);
  return eta_var.adj();
}

}  // namespace internal

/**
 * @tparam F
 * @tparam Theta
 * @tparam Eta
 * @tparam TupleArgs
 * @param f
 * @param theta
 * @param eta
 * @param ll_tup
 * @param msgs
 */
template <typename F, typename Theta, typename Eta, typename TupleArgs,
          require_eigen_vector_t<Theta>* = nullptr,
          require_eigen_t<Eta>* = nullptr,
          require_tuple_t<TupleArgs>* = nullptr>
inline auto log_likelihood(F&& f, const Theta& theta, const Eta& eta,
                           TupleArgs&& ll_tup, std::ostream* msgs) {
  return apply(
      [](auto&& f, auto&& theta, auto&& eta, auto&& msgs, auto&&... args) {
        return internal::log_likelihood(f, theta, eta, args..., msgs);
      },
      ll_tup, f, theta, eta, msgs);
}

/**
 * @tparam F
 * @tparam Theta
 * @tparam Eta
 * @tparam TupleArgs
 * @param f
 * @param theta
 * @param eta
 * @param gradient
 * @param hessian_block_size
 * @param ll_tuple
 * @param msgs
 */
template <typename F, typename Theta, typename Eta, typename TupleArgs,
          require_eigen_vector_t<Theta>* = nullptr,
          require_eigen_t<Eta>* = nullptr,
          require_tuple_t<TupleArgs>* = nullptr>
inline Eigen::SparseMatrix<double> diff(F&& f, const Theta& theta,
                                        const Eta& eta,
                                        Eigen::VectorXd& gradient,
                                        const Eigen::Index hessian_block_size,
                                        TupleArgs&& ll_tuple,
                                        std::ostream* msgs) {
  return apply(
      [](auto&& f, auto&& theta, auto&& eta, auto&& gradient,
         auto hessian_block_size, auto* msgs, auto&&... args) {
        return internal::diff(f, theta, eta, gradient, hessian_block_size,
                              args..., msgs);
      },
      ll_tuple, f, theta, eta, gradient, hessian_block_size, msgs);
}

/**
 * @tparam F
 * @tparam Theta
 * @tparam Eta
 * @tparam TupleArgs
 * @param f
 * @param theta
 * @param eta
 * @param ll_args
 * @param msgs
 */
template <typename F, typename Theta, typename Eta, typename TupleArgs,
          require_eigen_vector_t<Theta>* = nullptr,
          require_eigen_t<Eta>* = nullptr,
          require_tuple_t<TupleArgs>* = nullptr>
inline Eigen::VectorXd third_diff(F&& f, const Theta& theta, const Eta& eta,
                                  TupleArgs&& ll_args, std::ostream* msgs) {
  return apply(
      [](auto&& f, auto&& theta, auto&& eta, auto&& msgs, auto&&... args) {
        return internal::third_diff(f, theta, eta, args..., msgs);
      },
      ll_args, f, theta, eta, msgs);
}

/**
 * @tparam F
 * @tparam Theta
 * @tparam Eta
 * @tparam TupleArgs
 * @param f
 * @param theta
 * @param eta
 * @param A
 * @param hessian_block_size
 * @param ll_args
 * @param msgs
 */
template <typename F, typename Theta, typename Eta, typename TupleArgs,
          require_eigen_vector_t<Theta>* = nullptr,
          require_eigen_t<Eta>* = nullptr,
          require_tuple_t<TupleArgs>* = nullptr>
inline Eigen::VectorXd compute_s2(F&& f, const Theta& theta, const Eta& eta,
                                  const Eigen::MatrixXd& A,
                                  int hessian_block_size, TupleArgs&& ll_args,
                                  std::ostream* msgs) {
  return apply(
      [](auto&& f, auto&& theta, auto&& eta, auto&& A, auto hessian_block_size,
         auto* msgs, auto&&... args) {
        return internal::compute_s2(f, theta, eta, A, hessian_block_size,
                                    args..., msgs);
      },
      ll_args, f, theta, eta, A, hessian_block_size, msgs);
}

/**
 * @tparam F
 * @tparam V_t
 * @tparam Theta
 * @tparam Eta
 * @tparam TupleArgs
 * @param f
 * @param v
 * @param theta
 * @param eta
 * @param ll_args
 * @param msgs
 */
template <typename F, typename V_t, typename Theta, typename Eta,
          typename TupleArgs, require_tuple_t<TupleArgs>* = nullptr,
          require_eigen_vector_t<Theta>* = nullptr,
          require_eigen_t<Eta>* = nullptr>
inline plain_type_t<Eta> diff_eta_implicit(F&& f, const V_t& v,
                                           const Theta& theta, const Eta& eta,
                                           TupleArgs&& ll_args,
                                           std::ostream* msgs) {
  return apply(
      [](auto&& f, auto&& v, auto&& theta, auto&& eta, auto&& msgs,
         auto&&... args) {
        return internal::diff_eta_implicit(f, v, theta, eta, args..., msgs);
      },
      ll_args, f, v, theta, eta, msgs);
}

}  // namespace laplace_likelihood

}  // namespace math
}  // namespace stan

#endif
