#ifndef STAN_MATH_MIX_FUNCTOR_LAPLACE_LIKELIHOOD_HPP
#define STAN_MATH_MIX_FUNCTOR_LAPLACE_LIKELIHOOD_HPP

// #include <stan/math/mix/laplace/hessian_times_vector.hpp>
#include <stan/math/mix/functor/hessian_block_diag.hpp>
#include <stan/math/prim/fun.hpp>
#include <stan/math/prim/functor.hpp>
#include <Eigen/Sparse>

namespace stan {
namespace math {
  inline std::basic_ostream<char>* value_of(std::basic_ostream<char>*& pstream) { return pstream; }

/**
 * functions to compute the log density, first, second,
 * and third-order derivatives for a likelihoood specified by the user.
 */
namespace laplace_likelihood {
namespace internal {
/**
 * @tparam F Type of log likelihood function.
 * @tparam Theta Type of latent Gaussian variable.
 * @tparam Args Type of variadic arguments.
 * @param f Log likelihood function.
 * @param theta Latent Gaussian variable.
 * @param args Additional variational arguments for likelihood function.
 */
template <typename F, typename Theta, typename... Args,
          require_eigen_vector_t<Theta>* = nullptr>
inline auto log_likelihood(F&& f, Theta&& theta, Args&&... args) {
  return std::forward<F>(f)(std::forward<Theta>(theta),
                            std::forward<Args>(args)...);
}

template <typename PromotedType, typename... Args>
inline auto deep_copy_and_promote(Args&&... args) {
  return apply_if<is_any_var_scalar>([](auto&& arg){
    if constexpr (is_tuple<std::decay_t<decltype(arg)>>::value) {
      return deep_copy_and_promote<PromotedType>(std::forward<decltype(arg)>(arg));
    } else {
      return stan::math::eval(promote_scalar<PromotedType>(value_of(std::forward<decltype(arg)>(arg))));
    }
  }, args...);
}

template <typename EigVec, typename Tuple>
inline auto serialize_adjoints_impl(EigVec& vec, Tuple&& tup, std::size_t& i = 0) {
  for_each([&i, &vec](auto&& arg) {
    using arg_t = std::decay_t<decltype(arg)>;
    if constexpr (is_any_var_scalar<arg_t>::value) {
      if constexpr (is_stan_scalar<arg_t>::value) {
        vec(i++) = arg.adj();
        return;
      } else if constexpr (is_eigen<arg_t>::value || (is_std_vector<arg_t>::value &&
        !is_container<value_type_t<arg_t>>::value)) {
        for (int j = 0; j < arg.size(); ++j) {
          vec(i++) = arg(j).adj();
        }
        return;
      } else if constexpr (is_std_vector<arg_t>::value && is_container<value_type_t<arg_t>>::value) {
        for (int j = 0; j < arg.size(); ++j) {
          serialize_adjoints_impl(vec, std::forward_as_tuple(arg[j]), i);
        }
        return;
      } else if constexpr (is_tuple<arg_t>::value) {
        serialize_adjoints_impl(vec, arg, i);
        return;
      } else {
        static_assert(0, "Unsupported type");
      }
    }
  }, tup);
}
template <typename EigVec, typename Tuple>
inline auto serialize_adjoints(EigVec& vec, Tuple&& tup) {
  std::size_t i = 0;
  serialize_adjoints_impl(vec, tup, i);
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
  auto [theta_gradient, eta_gradient] = [&theta, &f](auto&&... args) {
    nested_rev_autodiff nested;
    Matrix<var, Dynamic, 1> theta_var = theta;
    auto hard_copy_args = deep_copy_and_promote<var>(args...);
    var f_var = stan::math::apply([](auto&& f, auto&& theta_var, auto&&... inner_args) {
      return f(theta_var, inner_args...);
    }, hard_copy_args, f, theta_var);
    grad(f_var.vi_);
    return std::make_pair(theta_var.adj().eval(),
      stan::math::filter<is_any_var_scalar>([](auto&& arg){
        return stan::math::eval(get_adj(std::forward<decltype(arg)>(arg)));
      }, hard_copy_args));
  }(args...);
  if (hessian_block_size == 1) {
    Eigen::VectorXd v = Eigen::VectorXd::Ones(theta_size);
    Eigen::VectorXd hessian_v = hessian_times_vector(f, theta, v, value_of(args)...);
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
        (-hessian_block_diag(f, theta, hessian_block_size, value_of(args)...))
            .eval());
  }
}

/**
 * @tparam F Type of log likelihood function.
 * @tparam Theta Type of latent Gaussian variable.
 * @tparam Args Type of variadic arguments for likelihood function.
 * @param f Log likelihood function.
 * @param theta Latent Gaussian variable.
 * @param args Variadic arguments for likelihood function.
 */
template <typename F, typename Theta, typename... Args,
          require_eigen_vector_t<Theta>* = nullptr>
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
 * @tparam Args Type of variadic arguments for likelihood function.
 * @param f Log likelihood function.
 * @param theta Latent Gaussian variable.
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
  Matrix<var, Dynamic, 1> theta_var = theta;
  int n_blocks = theta_size / hessian_block_size;
  fvar<fvar<var>> target_ffvar = 0;
  VectorXd v(theta_size);
  VectorXd w(theta_size);
  auto copy_vargs = deep_copy_and_promote<var>(args...);
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
    auto hard_copy_args = stan::math::apply_if<is_any_var_scalar>([](auto&& arg){
      return arg.template cast<fvar<var>>().eval();
    }, copy_vargs);
    target_ffvar += stan::math::apply([](auto&& f, auto&& theta_ffvar, auto&&... inner_args) {
      return f(theta_ffvar, inner_args...);
    }, hard_copy_args, f, theta_ffvar);
  }
  // TODO:(Steve) Use many small grad instead of the large one?
  grad(target_ffvar.d_.d_.vi_);
  auto eta_grad = stan::math::filter<is_any_var_scalar>([](auto&& arg){
      return stan::math::eval(0.5 * get_adj(std::forward<decltype(arg)>(arg)));
    }, copy_vargs);
  return std::make_pair((0.5 * theta_var.adj()).eval(),
    eta_grad);
}

/**
 * @tparam F Type of log likelihood function.
 * @tparam Theta Type of latent Gaussian variable.
 * @tparam Args Type of variational arguments for likelhood function.
 * @param f Log likelihood function.
 * @param v Initial tangent.
 * @param theta Latent Gaussian variable.
 * @param args Variadic arguments for likelhood function.
 */
template <typename F, typename V_t, typename Theta,
          typename... Args, require_eigen_vector_t<Theta>* = nullptr>
inline auto constexpr diff_eta_implicit(F&& f, const V_t& v,
                                           const Theta& theta,
                                           Args&&... args) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::VectorXd;
  constexpr bool contains_var = is_any_var_scalar<Args...>::value;
  if constexpr (!contains_var) {
    return std::make_tuple();
  }
  nested_rev_autodiff nested;
  auto copy_vargs = deep_copy_and_promote<var>(args...);
  // CHECK -- can we avoid declaring theta as fvar<var>?
  const Eigen::Index theta_size = theta.size();
  Matrix<var, Dynamic, 1> theta_var = theta;
  Matrix<fvar<var>, Dynamic, 1> theta_fvar(theta_size);
  for (Eigen::Index i = 0; i < theta_size; i++) {
    theta_fvar(i) = fvar<var>(theta_var(i), v(i));
  }

  auto hard_copy_args = stan::math::apply_if<is_any_var_scalar>([](auto&& arg){
      return arg.template cast<fvar<var>>();
  }, copy_vargs);
  fvar<var> f_fvar = stan::math::apply([](auto&& f, auto&& theta_fvar, auto&&... inner_args) {
    return f(theta_fvar, inner_args...);
  }, hard_copy_args, f, theta_fvar);
  grad(f_fvar.d_.vi_);
  return stan::math::filter<is_any_var_scalar>([](auto&& arg){
      return stan::math::eval(get_adj(std::forward<decltype(arg)>(arg)));
    }, copy_vargs);
}

}  // namespace internal

/**
 * @tparam F Type of log likelihood function.
 * @tparam Theta Type of latent Gaussian variable.
 * @tparam TupleArgs Type of arguments for covariance function.
 * @param f Log likelihood function.
 * @param theta Latent Gaussian model.
 * @param ll_tup Arguments for covariance function.
 * @param msgs stream messages.
 */
template <typename F, typename Theta, typename TupleArgs,
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
 * @tparam TupleArgs Type of arguments for covariance function.
 * @param f Log likelihood function.
 * @param theta Latent Gaussian model.
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
 * @tparam TupleArgs Type of arguments for covariance function.
 * @param f Log likelihood function.
 * @param theta Latent Gaussian variable.
 * @param ll_args Variadic arguments for likelihood function.
 * @param msgs Streaming message.
 */
template <typename F, typename Theta, typename TupleArgs,
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
 * @tparam TupleArgs Type of arguments for covariance function.
 * @param f Log likelihood function.
 * @param theta Latent Gaussian variable.
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
 * @tparam TupleArgs Type of variadic arguments for likelihood function.
 * @param f Log likelihood function.
 * @param v Initial tangent.
 * @param theta Latent Gaussian variable.
 * @param ll_args Variadic arguments for likelihood function.
 * @param msgs Streaming messages.
 */
template <typename F, typename V_t, typename Theta,
          typename TupleArgs, require_tuple_t<TupleArgs>* = nullptr,
          require_eigen_vector_t<Theta>* = nullptr>
inline auto diff_eta_implicit(F&& f, const V_t& v,
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
