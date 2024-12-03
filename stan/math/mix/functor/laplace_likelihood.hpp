#ifndef STAN_MATH_MIX_FUNCTOR_LAPLACE_LIKELIHOOD_HPP
#define STAN_MATH_MIX_FUNCTOR_LAPLACE_LIKELIHOOD_HPP

// #include <stan/math/mix/laplace/hessian_times_vector.hpp>
#include <stan/math/mix/functor/hessian_block_diag.hpp>
#include <stan/math/prim/fun.hpp>
#include <stan/math/prim/functor.hpp>
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
 * @tparam F Type of log likelihood function.
 * @tparam Theta Type of latent Gaussian variable.
 * @tparam Eta Type of parameter arguments for log likelihood function.
 * @tparam Args Type of variadic arguments.
 * @param f Log likelihood function.
 * @param theta Latent Gaussian variable.
 * @param eta Parameter arguments for likelihood function.
 * @param args Additional variational arguments for likelihood function.
 */
template <typename F, typename Theta, typename Eta, typename... Args,
          require_eigen_vector_t<Theta>* = nullptr,
          require_eigen_t<Eta>* = nullptr>
inline auto log_likelihood(F&& f, Theta&& theta, Eta&& eta, Args&&... args) {
  return std::forward<F>(f)(std::forward<Theta>(theta), std::forward<Eta>(eta),
                            std::forward<Args>(args)...);
}

template <typename Tuple, require_tuple_t<Tuple>* = nullptr>
inline auto get_vector_adjoints(Tuple&& tup) {
  if constexpr (!is_any_arg_ad_v<Tuple>) {
    return Eigen::Matrix<double, 0, 0>{};
  } else {
    const Eigen::Index num_vars = count_vars(tup);
    Eigen::Matrix<double, -1, 1> adjoints = Eigen::Matrix<double, -1, 1>::Zeros(num_vars);
    accumulate_adjoints(adjoints.data(), tup);
  }
}

template <typename... Args>
inline auto copy_vars(Args&&... args) {
  return filter_map<is_scalar_var>([](auto&& arg) {
      // TODO:(Steve) Needs to handle std::vector
      return std::decay_t<decltype(arg)>(value_of(arg));
    }, std::forward<Args>(args)...);;
}

/**
 * @tparam F Type of log likelihood function.
 * @tparam Theta Type of latent Gaussian variable.
 * @tparam Args Type of variadic arguments.
 * @param f Log likelihood function.
 * @param theta Latent Gaussian model.
 * @param gradient Gradient of likelihood returned by function.
 * @param hessian_block_size If the Hessian of the log likelihood function w.r.t
 *                           the latent Gaussian variable is block-diagonal,
 *                           size of each block.
 * @param args Variadic arguments for the likelihood function.
 */
template <typename F, typename Theta, typename... Args,
          require_eigen_vector_t<Theta>* = nullptr>
inline auto diff(F&& f, const Theta& theta,
                 const Eigen::Index hessian_block_size, Args&&... args) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  const Eigen::Index theta_size = theta.size();
  Eigen::Matrix<double, Theta::RowsAtCompileTime, Theta::ColsAtCompileTime>
      theta_gradient;
  const auto num_vars = count_vars(args...);
  // TODO:(Steve) Write is any_autodiff_v
  constexpr bool is_any_arg_ad = is_any_autodiff_v<Args...>;
  Eigen::Matrix<double, is_any_arg_ad ? -1 : 0 , is_any_arg_ad ? 1 : 0>
      eta_gradient(num_vars);
  {
    nested_rev_autodiff nested;
    Matrix<var, Dynamic, 1> theta_var = theta;
    auto var_copied_args = copy_vars(std::forward<Args>(args)...);
    var f_var = stan::math::apply([](f&& f, auto&& theta_var, auto&&... inner_args) {
      return f(theta_var, inner_args...);
    }, var_copied_args, f, theta_var);
    grad(f_var.vi_);
    theta_gradient = theta_var.adj();
    get_adjoints(eta_gradient, var_copied_args);
  }

  if (hessian_block_size == 1) {
    Eigen::VectorXd v = Eigen::VectorXd::Ones(theta_size);
    Eigen::VectorXd hessian_v = hessian_times_vector(f, theta, v, args...);
    Eigen::SparseMatrix<double> hessian_theta(theta_size, theta_size);
    hessian_theta.reserve(Eigen::VectorXi::Constant(theta_size, 1));
    for (Eigen::Index i = 0; i < theta_size; i++) {
      hessian_theta.insert(i, i) = hessian_v(i);
    }
    return std::make_tuple(std::move(theta_gradient), std::move(eta_gradient),
                           (-hessian_theta).eval());
  } else {
    return std::make_tuple(
        std::move(theta_gradient), std::move(eta_gradient),
        (-hessian_block_diag(f, theta, hessian_block_size, args...))
            .eval());
  }
}

/**
 * @tparam F Type of log likelihood function.
 * @tparam Theta Type of latent Gaussian variable.
 * @tparam Eta Type of parameter arguments for likelihood function.
 * @tparam Args Type of variadic arguments for likelihood function.
 * @param f Log likelihood function.
 * @param theta Latent Gaussian variable.
 * @param eta Parameter arguments for likelihood function.
 * @param args Variadic arguments for likelihood function.
 */
template <typename F, typename Theta, typename... Args,
          require_eigen_vector_t<Theta>* = nullptr,
          require_t<is_all_arithmetic_scalar<Args>>* = nullptr>
inline Eigen::VectorXd third_diff(F&& f, const Theta& theta,
                                  Args&&... args) {
  nested_rev_autodiff nested;
  const Eigen::Index theta_size = theta.size();
  Eigen::Matrix<var, Eigen::Dynamic, 1> theta_var = theta;
  Eigen::Matrix<fvar<fvar<var>>, Eigen::Dynamic, 1> theta_ffvar(theta_size);
  for (Eigen::Index i = 0; i < theta_size; ++i) {
    theta_ffvar(i) = fvar<fvar<var>>(fvar<var>(theta_var(i), 1.0), 1.0);
  }
  fvar<fvar<var>> ftheta_ffvar = f(theta_ffvar, args...);
  grad(ftheta_ffvar.d_.d_.vi_);
  return theta_var.adj();
}

/**
 * @tparam F Type of log likelihood function.
 * @tparam Theta Type of latent Gaussian variable.
 * @tparam Eta Type of parameter arguments for likelihood function.
 * @tparam Args Type of variadic arguments for likelihood function.
 * @param f Log likelihood function.
 * @param theta Latent Gaussian variable.
 * @param eta Parameter arguments for
 * @param A Matrix storing initial tangents for higher-order differentiation
 *        (line 21 in Algorithm 4, https://arxiv.org/pdf/2306.14976)
 * @param hessian_block_size If the Hessian of the log likelihood w.r.t theta
 *                           is block diagonal, size of each block.
 * @param args Variational arguments for likelihood function.
 */
template <typename F, typename Theta, typename... Args,
          require_eigen_vector_t<Theta>* = nullptr>
inline auto compute_s2(F&& f, const Theta& theta,
                       const Eigen::MatrixXd& A, const int hessian_block_size,
                       Args&&... args) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::MatrixXd;
  using Eigen::VectorXd;

  nested_rev_autodiff nested;
  const Eigen::Index theta_size = theta.size();
  const Eigen::Index parm_size = theta_size + eta_size;
  Matrix<var, Dynamic, 1> theta_var = theta;
  const auto num_vars = count_vars(args...);
  // TODO:(Steve) Write is any_autodiff_v
  constexpr bool is_any_arg_ad = is_any_autodiff_v<Args...>;
  Eigen::Matrix<double, is_any_arg_ad ? -1 : 0 , is_any_arg_ad ? 1 : 0>
      eta_gradient(num_vars);
  auto var_copied_args = copy_vars(std::forward<Args>(args)...);

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
    Matrix<fvar<fvar<var>>, Dynamic, 1> theta_ffvar(theta_size);
    for (int j = 0; j < theta_size; ++j) {
      theta_ffvar(j) = fvar<fvar<var>>(theta_fvar(j), w(j));
    }
    auto ffvar_copied_args = stan::math::filter_map<is_scalar_var>([](auto&& arg) {
      return promote_scalar_t<fvar<fvar<var>>>(arg);
    }, var_copied_args);
    target_ffvar += stan::math::apply([](auto&& f, auto&& theta_ffvar, auto&&... inner_ffvars) {
      return f(theta_ffvar, inner_ffvars...);
    }, ffvar_copied_args, f, theta_ffvar);
  }
  grad(target_ffvar.d_.d_.vi_);
  get_adjoints(eta_gradient, var_copied_args);
  return std::make_pair((0.5 * theta_var.adj()).eval(),
                        (0.5 * eta_gradient).eval());
}

/**
 * @tparam F Type of log likelihood function.
 * @tparam Theta Type of latent Gaussian variable.
 * @tparam Eta Type of parameter arguments for likelhood function.
 * @tparam Args Type of variational arguments for likelhood function.
 * @param f Log likelihood function.
 * @param v Initial tangent.
 * @param theta Latent Gaussian variable.
 * @param eta Parameter arguments for likelhood function.
 * @param args Variadic arguments for likelhood function.
 */
template <typename F, typename V_t, typename Theta,
          typename... Args, require_eigen_vector_t<Theta>* = nullptr,
          require_eigen_t<Eta>* = nullptr>
inline plain_type_t<Eta> diff_eta_implicit(F&& f, const V_t& v,
                                           const Theta& theta,
                                           Args&&... args) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::VectorXd;
  if constexpr (Eta::RowsAtCompileTime == 0 && Eta::ColsAtCompileTime == 0) {
    return Eigen::Matrix<double, 0, 0>{};
  }
  nested_rev_autodiff nested;
  const auto num_vars = count_vars(args...);
  // TODO:(Steve) Write is any_autodiff_v
  constexpr bool is_any_arg_ad = is_any_autodiff_v<Args...>;
  Eigen::Matrix<double, is_any_arg_ad ? -1 : 0 , is_any_arg_ad ? 1 : 0>
      eta_gradient(num_vars);
  auto var_copied_args = copy_vars(std::forward<Args>(args)...);
  // CHECK -- can we avoid declaring theta as fvar<var>?
  // We currently compute derivatives wrt eta, which is not needed.
  const Eigen::Index theta_size = theta.size();
  Matrix<var, Dynamic, 1> theta_var = theta;
  Matrix<fvar<var>, Dynamic, 1> theta_fvar(theta_size);
  for (Eigen::Index i = 0; i < theta_size; i++) {
    theta_fvar(i) = fvar<var>(theta_var(i), v(i));
  }
  auto fvar_copied_args = stan::math::filter_map<is_scalar_var>([](auto&& arg) {
    return promote_scalar_t<fvar<var>>(arg);
  }, var_copied_args);

  fvar<var> f_fvar = stan::math::apply([](auto&& f, auto&& theta_fvar, auto&&... inner_fvars) {
    return f(theta_ffvar, inner_fvars...);
  }, fvar_copied_args, f, theta_fvar);
  grad(f_fvar.d_.vi_);
  get_adjoints(eta_gradient, var_copied_args);
  return eta_gradient;
}

}  // namespace internal

/**
 * @tparam F Type of log likelihood function.
 * @tparam Theta Type of latent Gaussian variable.
 * @tparam Eta Type of parameter argument for likelihood function.
 * @tparam TupleArgs Type of arguments for covariance function.
 * @param f Log likelihood function.
 * @param theta Latent Gaussian model.
 * @param eta Parameter arguments for likelihood function.
 * @param ll_tup Arguments for covariance function.
 * @param msgs stream messages.
 */
template <typename F, typename Theta, typename Eta, typename TupleArgs,
          require_eigen_vector_t<Theta>* = nullptr,
          require_tuple_t<TupleArgs>* = nullptr>
inline auto log_likelihood(F&& f, const Theta& theta,
                           TupleArgs&& ll_tup, std::ostream* msgs) {
  return apply(
      [](auto&& f, auto&& theta, auto&& msgs, auto&&... args) {
        return internal::log_likelihood(f, theta, args..., msgs);
      },
      ll_tup, f, theta, msgs);
}

/**
 * @tparam F Type of log likelihood function.
 * @tparam Theta Type of latent Gaussian variable.
 * @tparam Eta Type of parameter variable.
 * @tparam TupleArgs Type of arguments for covariance function.
 * @param f Log likelihood function.
 * @param theta Latent Gaussian model.
 * @param eta Parameter arguments for likelihood function.
 * @param gradient Vector to store gradient of log likelihood w.r.t theta.
 * @param hessian_block_size If Hessian of log likelihood w.r.t theta is
 *                           block diagonal, size of block.
 * @param ll_tuple Arguments of covariance function.
 * @param msgs Stream messages.
 */
template <typename F, typename Theta, typename TupleArgs,
          require_eigen_vector_t<Theta>* = nullptr,
          require_tuple_t<TupleArgs>* = nullptr>
inline auto diff(F&& f, const Theta& theta,
                 const Eigen::Index hessian_block_size, TupleArgs&& ll_tuple,
                 std::ostream* msgs) {
  return apply(
      [](auto&& f, auto&& theta, auto hessian_block_size,
         auto* msgs, auto&&... args) {
        return internal::diff(f, theta, hessian_block_size, args..., msgs);
      },
      ll_tuple, f, theta, hessian_block_size, msgs);
}

/**
 * @tparam F Type of log likelhood function.
 * @tparam Theta Type of latent Gaussian variable.
 * @tparam Eta Type of parameter arguments for likelhood function.
 * @tparam TupleArgs Type of arguments for covariance function.
 * @param f Log likelihood function.
 * @param theta Latent Gaussian variable.
 * @param eta Parameter argument for likelihood funciton.
 * @param ll_args Variadic arguments for likelihood function.
 * @param msgs Streaming message.
 */
template <typename F, typename Theta, typename Eta, typename TupleArgs,
          require_eigen_vector_t<Theta>* = nullptr,
          require_tuple_t<TupleArgs>* = nullptr>
inline Eigen::VectorXd third_diff(F&& f, const Theta& theta,
                                  TupleArgs&& ll_args, std::ostream* msgs) {
  return apply(
      [](auto&& f, auto&& theta, auto&& msgs, auto&&... args) {
        return internal::third_diff(f, theta, args..., msgs);
      },
      ll_args, f, theta, msgs);
}

/**
 * @tparam F Type of log likelhood function.
 * @tparam Theta Type of latent Gaussian ba
 * @tparam Eta Type of parameter argument for likelihood function.
 * @tparam TupleArgs Type of arguments for covariance function.
 * @param f Log likelihood function.
 * @param theta Latent Gaussian variable.
 * @param eta Parameter  arguments for likelihood function.
 * @param A Matrix storing initial tangents for higher-order differentiation
 *        (line 21 in Algorithm 4, https://arxiv.org/pdf/2306.14976)
 * @param hessian_block_size If Hessian of log likelihood w.r.t theta is
 *                           block diagonal, size of block.
 * @param ll_args Variadic arguments for likelihood function.
 * @param msgs Streaming messages.
 */
template <typename F, typename Theta, typename TupleArgs,
          require_eigen_vector_t<Theta>* = nullptr,
          require_tuple_t<TupleArgs>* = nullptr>
inline auto compute_s2(F&& f, const Theta& theta,
                       const Eigen::MatrixXd& A, int hessian_block_size,
                       TupleArgs&& ll_args, std::ostream* msgs) {
  return apply(
      [](auto&& f, auto&& theta, auto&& A, auto hessian_block_size,
         auto* msgs, auto&&... args) {
        return internal::compute_s2(f, theta, A, hessian_block_size,
                                    args..., msgs);
      },
      ll_args, f, theta, A, hessian_block_size, msgs);
}

/**
 * @tparam F Type of log likelihood function.
 * @tparam V_t Type of initial tangent.
 * @tparam Theta Type of latent Gaussian variable.
 * @tparam Eta Type of parameter arguments for likelihood function.
 * @tparam TupleArgs Type of variadic arguments for likelihood function.
 * @param f Log likelihood function.
 * @param v Initial tangent.
 * @param theta Latent Gaussian variable.
 * @param eta Parameter argument for likelihood function.
 * @param ll_args Variadic arguments for likelihood function.
 * @param msgs Streaming messages.
 */
template <typename F, typename V_t, typename Theta,
          typename TupleArgs, require_tuple_t<TupleArgs>* = nullptr,
          require_eigen_vector_t<Theta>* = nullptr>
inline plain_type_t<Eta> diff_eta_implicit(F&& f, const V_t& v,
                                           const Theta& theta,
                                           TupleArgs&& ll_args,
                                           std::ostream* msgs) {
  return apply(
      [](auto&& f, auto&& v, auto&& theta, auto&& msgs,
         auto&&... args) {
        return internal::diff_eta_implicit(f, v, theta, args..., msgs);
      },
      ll_args, f, v, theta, msgs);
}

}  // namespace laplace_likelihood

}  // namespace math
}  // namespace stan

#endif
