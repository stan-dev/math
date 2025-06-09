#ifndef STAN_MATH_MIX_FUNCTOR_LAPLACE_LIKELIHOOD_HPP
#define STAN_MATH_MIX_FUNCTOR_LAPLACE_LIKELIHOOD_HPP

#include <stan/math/mix/functor/hessian_block_diag.hpp>
#include <stan/math/prim/functor.hpp>
#include <stan/math/prim/fun.hpp>

namespace stan {
namespace math {

/**
 * functions to compute the log density, first, second,
 * and third-order derivatives for a likelihoood specified by the user.
 */
namespace laplace_likelihood {
namespace internal {
/**
 * @tparam F A functor with `opertor()(Args&&...)` returning a scalar
 * @tparam Theta A class assignable to an Eigen vector type
 * @tparam Stream Type of stream for messages.
 * @tparam Args Type of variadic arguments.
 * @param f Log likelihood function.
 * @param theta Latent Gaussian variable.
 * @param msgs Stream for messages.
 * @param args Additional variational arguments for likelihood function.
 */
template <typename F, typename Theta, typename Stream, typename... Args,
          require_eigen_vector_t<Theta>* = nullptr>
inline auto log_likelihood(F&& f, Theta&& theta, Stream* msgs, Args&&... args) {
  return std::forward<F>(f)(std::forward<Theta>(theta),
                            std::forward<Args>(args)..., msgs);
}

/**
 * Decide if object should be deep or shallow copied when
 * using @ref conditional_copy_and_promote .
 */
enum class COPY_TYPE { SHALLOW = 0, DEEP = 1 };

/**
 * Conditional copy and promote a type's scalar type to a `PromotedType`.
 * @tparam Filter type trait with a static constexpr bool member `value`
 *  that is true if the type should be promoted. Otherwise, the type is
 *  left unchanged.
 * @tparam PromotedType type to promote the scalar to.
 * @tparam CopyType type of copy to perform.
 * @tparam Args variadic arguments.
 * @param args variadic arguments to conditionally copy and promote.
 * @return a tuple where each element is either a reference to the original
 * argument or a promoted copy of the argument.
 */
template <template <typename...> class Filter,
          typename PromotedType = stan::math::var,
          COPY_TYPE CopyType = COPY_TYPE::DEEP, typename... Args>
inline auto conditional_copy_and_promote(Args&&... args) {
  return map_if<Filter>(
      [](auto&& arg) {
        if constexpr (is_tuple_v<decltype(arg)>) {
          return stan::math::apply(
              [](auto&&... inner_args) {
                return make_holder_tuple(
                    conditional_copy_and_promote<Filter, PromotedType,
                                                 CopyType>(
                        std::forward<decltype(inner_args)>(inner_args))...);
              },
              std::forward<decltype(arg)>(arg));
        } else if constexpr (is_std_vector_v<decltype(arg)>) {
          std::vector<decltype(conditional_copy_and_promote<
                               Filter, PromotedType, CopyType>(arg[0]))>
              ret;
          for (std::size_t i = 0; i < arg.size(); ++i) {
            ret.push_back(
                conditional_copy_and_promote<Filter, PromotedType, CopyType>(
                    arg[i]));
          }
          return ret;
        } else {
          if constexpr (CopyType == COPY_TYPE::DEEP) {
            return stan::math::eval(promote_scalar<PromotedType>(
                value_of_rec(std::forward<decltype(arg)>(arg))));
          } else if (CopyType == COPY_TYPE::SHALLOW) {
            if constexpr (std::is_same_v<PromotedType,
                                         scalar_type_t<decltype(arg)>>) {
              return std::forward<decltype(arg)>(arg);
            } else {
              return stan::math::eval(promote_scalar<PromotedType>(
                  std::forward<decltype(arg)>(arg)));
            }
          }
        }
      },
      std::forward<Args>(args)...);
}

template <typename PromotedType, typename... Args>
inline auto deep_copy_vargs(Args&&... args) {
  return conditional_copy_and_promote<is_any_var_scalar, PromotedType,
                                      COPY_TYPE::DEEP>(
      std::forward<Args>(args)...);
}

template <typename PromotedType, typename... Args>
inline auto shallow_copy_vargs(Args&&... args) {
  return conditional_copy_and_promote<is_any_var_scalar, PromotedType,
                                      COPY_TYPE::SHALLOW>(
      std::forward<Args>(args)...);
}

/**
 * Computes theta gradient and negative block diagonal Hessian of `f` wrt
 * `theta` and `args...`
 * @note If `Args` contains \ref var types then their adjoints will be
 * calculated as a side effect.
 * @tparam F A functor with `opertor()(Args&&...)` returning a scalar
 * @tparam Theta A class assignable to an Eigen vector type
 * @tparam Stream Type of stream for messages.
 * @tparam Args Type of variadic arguments.
 * @param f Log likelihood function.
 * @param theta Latent Gaussian model.
 * @param hessian_block_size If the Hessian of the log likelihood function w.r.t
 *                           the latent Gaussian variable is block-diagonal,
 *                           size of each block.
 * @param msgs Stream for messages.
 * @param args Variadic arguments for the likelihood function.
 */
template <typename F, typename Theta, typename Stream, typename... Args,
          require_eigen_vector_vt<std::is_arithmetic, Theta>* = nullptr>
inline auto diff(F&& f, Theta&& theta, const Eigen::Index hessian_block_size,
                 Stream* msgs, Args&&... args) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  const Eigen::Index theta_size = theta.size();
  auto theta_gradient = [&theta, &f, &msgs](auto&&... args) {
    nested_rev_autodiff nested;
    Matrix<var, Dynamic, 1> theta_var = theta;
    var f_var = f(theta_var, args..., msgs);
    grad(f_var.vi_);
    return theta_var.adj().eval();
  }(args...);
  if (hessian_block_size == 1) {
    auto v = Eigen::VectorXd::Ones(theta_size);
    Eigen::VectorXd hessian_v = Eigen::VectorXd::Zero(theta_size);
    hessian_times_vector(f, hessian_v, std::forward<Theta>(theta), std::move(v),
                         value_of(args)..., msgs);
    Eigen::SparseMatrix<double> hessian_theta(theta_size, theta_size);
    hessian_theta.reserve(Eigen::VectorXi::Constant(theta_size, 1));
    for (Eigen::Index i = 0; i < theta_size; i++) {
      hessian_theta.insert(i, i) = hessian_v(i);
    }
    return std::make_pair(std::move(theta_gradient), (-hessian_theta).eval());
  } else {
    return std::make_pair(
        std::move(theta_gradient),
        (-hessian_block_diag(f, std::forward<Theta>(theta), hessian_block_size,
                             value_of(args)..., msgs))
            .eval());
  }
}

/**
 * Compute third order derivative of `f` wrt `theta` and `args...`
 * @note If `Args` contains \ref var types then their adjoints will be
 * calculated as a side effect.
 * @tparam F A functor with `opertor()(Args&&...)` returning a scalar
 * @tparam Theta A class assignable to an Eigen vector type
 * @tparam Stream Type of stream for messages.
 * @tparam Args Type of variadic arguments for likelihood function.
 * @param f Log likelihood function.
 * @param theta Latent Gaussian variable.
 * @param msgs Stream for messages.
 * @param args Variadic arguments for likelihood function.
 */
template <typename F, typename Theta, typename Stream, typename... Args,
          require_eigen_vector_t<Theta>* = nullptr>
inline Eigen::VectorXd third_diff(F&& f, Theta&& theta, Stream&& msgs,
                                  Args&&... args) {
  nested_rev_autodiff nested;
  const Eigen::Index theta_size = theta.size();
  Eigen::Matrix<var, Eigen::Dynamic, 1> theta_var = std::forward<Theta>(theta);
  Eigen::Matrix<fvar<fvar<var>>, Eigen::Dynamic, 1> theta_ffvar(theta_size);
  for (Eigen::Index i = 0; i < theta_size; ++i) {
    theta_ffvar(i) = fvar<fvar<var>>(fvar<var>(theta_var(i), 1.0), 1.0);
  }
  fvar<fvar<var>> ftheta_ffvar = f(theta_ffvar, args..., msgs);
  grad(ftheta_ffvar.d_.d_.vi_);
  return theta_var.adj().eval();
}

/**
 * The derivative of the log likelihood wrt `theta` evaluated at the mode.
 * @brief Compute $s_2 = \Delta_{\theta} log \pi_G(y|\phi,\eta) = -\frac{1}{2}
 * trace((K^{-1}+W)^{-1})$
 * @note Equation 15 in https://arxiv.org/pdf/2306.14976
 * @note If `Args` contains \ref var types then their adjoints will be
 * calculated as a side effect.
 * @tparam F A functor with `opertor()(Args&&...)` returning a scalar
 * @tparam Theta An Eigen Matrix
 * @tparam AMat An Eigen Matrix
 * @tparam Stream Type of stream for messages.
 * @tparam Args Type of variadic arguments for likelihood function.
 * @param f Log likelihood function.
 * @param theta Latent Gaussian variable.
 * @param A Matrix storing initial tangents for higher-order differentiation
 *        (line 21 in Algorithm 4, https://arxiv.org/pdf/2306.14976)
 * @param hessian_block_size If the Hessian of the log likelihood w.r.t theta
 *                           is block diagonal, size of each block.
 * @param msgs Stream for messages.
 * @param args Variational arguments for likelihood function.
 */
template <typename F, typename Theta, typename AMat, typename Stream,
          typename... Args, require_eigen_vector_t<Theta>* = nullptr>
inline auto compute_s2(F&& f, Theta&& theta, AMat&& A,
                       const int hessian_block_size, Stream* msgs,
                       Args&&... args) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::MatrixXd;
  using Eigen::VectorXd;

  nested_rev_autodiff nested;
  const Eigen::Index theta_size = theta.size();
  Matrix<var, Dynamic, 1> theta_var = std::forward<Theta>(theta);
  int n_blocks = theta_size / hessian_block_size;
  VectorXd v(theta_size);
  VectorXd w(theta_size);
  Matrix<fvar<fvar<var>>, Dynamic, 1> theta_ffvar(theta_size);
  auto shallow_copy_args
      = shallow_copy_vargs<fvar<fvar<var>>>(std::forward_as_tuple(args...));
  for (Eigen::Index i = 0; i < hessian_block_size; ++i) {
    nested_rev_autodiff nested;
    v.setZero();
    for (int j = i; j < theta_size; j += hessian_block_size) {
      v(j) = 1;
    }
    w.setZero();
    for (int j = 0; j < n_blocks; ++j) {
      for (int k = 0; k < hessian_block_size; ++k) {
        w(k + j * hessian_block_size)
            = A(k + j * hessian_block_size, i + j * hessian_block_size);
      }
    }
    for (int j = 0; j < theta_size; ++j) {
      theta_ffvar(j) = fvar<fvar<var>>(fvar<var>(theta_var(j), v(j)), w(j));
    }
    fvar<fvar<var>> target_ffvar = stan::math::apply(
        [](auto&& f, auto&& theta_ffvar, auto&& msgs, auto&&... inner_args) {
          return f(theta_ffvar, inner_args..., msgs);
        },
        shallow_copy_args, f, theta_ffvar, msgs);
    grad(target_ffvar.d_.d_.vi_);
  }
  return (0.5 * theta_var.adj()).eval();
}

/**
 * Compute second order gradient of `f` wrt `theta` and `args...`
 * @note See proposition 2 in https://arxiv.org/pdf/2306.14976
 * See lines 31-37 in Algorithm 4
 * If `Args` contains \ref var types then their adjoints will be
 * calculated as a side effect.
 * @tparam F A functor with `opertor()(Args&&...)` returning a scalar
 * @tparam V_t A type assignable to an Eigen vector type
 * @tparam Theta A type assignable to an Eigen vector type
 * @tparam Stream Type of stream for messages.
 * @tparam Args Parameter pack of arguments to `F`'s `operator()`
 * @param f Log likelihood function.
 * @param v Initial tangent.
 * @param theta Latent Gaussian variable.
 * @param msgs Stream for messages.
 * @param args Variadic arguments for likelhood function.
 * @return `args` which are var types will have their adjoints set as a side
 * effect of this function.
 */
template <typename F, typename V_t, typename Theta, typename Stream,
          typename... Args, require_eigen_vector_t<Theta>* = nullptr>
inline auto diff_eta_implicit(F&& f, V_t&& v, Theta&& theta, Stream* msgs,
                              Args&&... args) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::VectorXd;
  constexpr bool contains_var = is_any_var_scalar<Args...>::value;
  if constexpr (!contains_var) {
    return;
  }
  nested_rev_autodiff nested;
  const Eigen::Index theta_size = theta.size();
  Matrix<var, Dynamic, 1> theta_var = std::forward<Theta>(theta);
  Matrix<fvar<var>, Dynamic, 1> theta_fvar(theta_size);
  for (Eigen::Index i = 0; i < theta_size; i++) {
    theta_fvar(i) = fvar<var>(theta_var(i), v(i));
  }
  auto shallow_copy_args
      = shallow_copy_vargs<fvar<var>>(std::forward_as_tuple(args...));
  fvar<var> f_sum = stan::math::apply(
      [](auto&& f, auto&& theta_fvar, auto&& msgs, auto&&... inner_args) {
        return f(theta_fvar, inner_args..., msgs);
      },
      shallow_copy_args, f, theta_fvar, msgs);
  grad(f_sum.d_.vi_);
}

}  // namespace internal

/**
 * A wrapper that accepts a tuple as arguments.
 * @tparam F A functor with `opertor()(Args&&...)` returning a scalar
 * @tparam Theta A class assignable to an Eigen vector type
 * @tparam TupleArgs Type of arguments for covariance function.
 * @tparam Stream Type of stream for messages.
 * @param f Log likelihood function.
 * @param theta Latent Gaussian model.
 * @param ll_tup Arguments for likelihood function
 * @param msgs stream messages.
 */
template <typename F, typename Theta, typename TupleArgs, typename Stream,
          require_eigen_vector_t<Theta>* = nullptr,
          require_tuple_t<TupleArgs>* = nullptr>
inline auto log_likelihood(F&& f, Theta&& theta, TupleArgs&& ll_tup,
                           Stream* msgs) {
  return apply(
      [](auto&& f, auto&& theta, auto&& msgs, auto&&... args) {
        return internal::log_likelihood(
            std::forward<decltype(f)>(f), std::forward<decltype(theta)>(theta),
            msgs, std::forward<decltype(args)>(args)...);
      },
      std::forward<TupleArgs>(ll_tup), std::forward<F>(f),
      std::forward<Theta>(theta), msgs);
}

/**
 * A wrapper that accepts a tuple as arguments.
 * @tparam F A functor with `opertor()(Args&&...)` returning a scalar
 * @tparam Theta A class assignable to an Eigen vector type
 * @tparam TupleArgs Type of arguments for covariance function.
 * @tparam Stream Type of stream for messages.
 * @param f Log likelihood function.
 * @param theta Latent Gaussian model.
 * @param hessian_block_size If Hessian of log likelihood w.r.t theta is
 *                           block diagonal, size of block.
 * @param ll_tuple Arguments for likelihood function
 * @param msgs Stream messages.
 */
template <typename F, typename Theta, typename TupleArgs, typename Stream,
          require_eigen_vector_t<Theta>* = nullptr,
          require_tuple_t<TupleArgs>* = nullptr>
inline auto diff(F&& f, Theta&& theta, const Eigen::Index hessian_block_size,
                 TupleArgs&& ll_tuple, Stream* msgs) {
  return apply(
      [](auto&& f, auto&& theta, auto hessian_block_size, auto* msgs,
         auto&&... args) {
        return internal::diff(
            std::forward<decltype(f)>(f), std::forward<decltype(theta)>(theta),
            hessian_block_size, msgs, std::forward<decltype(args)>(args)...);
      },
      std::forward<TupleArgs>(ll_tuple), std::forward<F>(f),
      std::forward<Theta>(theta), hessian_block_size, msgs);
}

/**
 * A wrapper that accepts a tuple as arguments.
 * @tparam F Type of log likelhood function.
 * @tparam Theta A class assignable to an Eigen vector type
 * @tparam TupleArgs Type of arguments for covariance function.
 * @tparam Stream Type of stream for messages.
 * @param f Log likelihood function.
 * @param theta Latent Gaussian variable.
 * @param ll_args Variadic arguments for likelihood function.
 * @param msgs Streaming message.
 */
template <typename F, typename Theta, typename TupleArgs, typename Stream,
          require_eigen_vector_t<Theta>* = nullptr,
          require_tuple_t<TupleArgs>* = nullptr>
inline Eigen::VectorXd third_diff(F&& f, Theta&& theta, TupleArgs&& ll_args,
                                  Stream* msgs) {
  return apply(
      [](auto&& f, auto&& theta, auto&& msgs, auto&&... args) {
        return internal::third_diff(std::forward<decltype(f)>(f),
                                    std::forward<decltype(theta)>(theta), msgs,
                                    std::forward<decltype(args)>(args)...);
      },
      std::forward<TupleArgs>(ll_args), std::forward<F>(f),
      std::forward<Theta>(theta), msgs);
}

/**
 * A wrapper that accepts a tuple as arguments.
 * @tparam F Type of log likelhood function.
 * @tparam Theta Type of latent Gaussian ba
 * @tparam TupleArgs Type of arguments for covariance function.
 * @tparam Stream Type of stream for messages.
 * @param f Log likelihood function.
 * @param theta Latent Gaussian variable.
 * @param A Matrix storing initial tangents for higher-order differentiation
 *        (line 21 in Algorithm 4, https://arxiv.org/pdf/2306.14976)
 * @param hessian_block_size If Hessian of log likelihood w.r.t theta is
 *                           block diagonal, size of block.
 * @param ll_args Variadic arguments for likelihood function.
 * @param msgs Streaming messages.
 */
template <typename F, typename Theta, typename AMat, typename TupleArgs,
          typename Stream, require_eigen_vector_t<Theta>* = nullptr,
          require_tuple_t<TupleArgs>* = nullptr>
inline auto compute_s2(F&& f, Theta&& theta, AMat&& A, int hessian_block_size,
                       TupleArgs&& ll_args, Stream* msgs) {
  return apply(
      [](auto&& f, auto&& theta, auto&& A, auto hessian_block_size, auto* msgs,
         auto&&... args) {
        return internal::compute_s2(
            std::forward<decltype(f)>(f), std::forward<decltype(theta)>(theta),
            std::forward<decltype(A)>(A), hessian_block_size, msgs,
            std::forward<decltype(args)>(args)...);
      },
      std::forward<TupleArgs>(ll_args), std::forward<F>(f),
      std::forward<Theta>(theta), std::forward<AMat>(A), hessian_block_size,
      msgs);
}

/**
 * A wrapper that accepts a tuple as arguments.
 * @tparam F A functor with `opertor()(Args&&...)` returning a scalar
 * @tparam V_t Type of initial tangent.
 * @tparam Theta A class assignable to an Eigen vector type
 * @tparam TupleArgs Type of variadic arguments for likelihood function.
 * @tparam Stream Type of stream for messages.
 * @param f Log likelihood function.
 * @param v Initial tangent.
 * @param theta Latent Gaussian variable.
 * @param ll_args Variadic arguments for likelihood function.
 * @param msgs Streaming messages.
 */
template <typename F, typename V_t, typename Theta, typename TupleArgs,
          typename Stream, require_tuple_t<TupleArgs>* = nullptr,
          require_eigen_vector_t<Theta>* = nullptr>
inline auto diff_eta_implicit(F&& f, V_t&& v, Theta&& theta,
                              TupleArgs&& ll_args, Stream* msgs) {
  return apply(
      [](auto&& f, auto&& v, auto&& theta, auto&& msgs, auto&&... args) {
        return internal::diff_eta_implicit(
            std::forward<decltype(f)>(f), std::forward<decltype(v)>(v),
            std::forward<decltype(theta)>(theta), msgs,
            std::forward<decltype(args)>(args)...);
      },
      std::forward<TupleArgs>(ll_args), std::forward<F>(f),
      std::forward<V_t>(v), std::forward<Theta>(theta), msgs);
}

}  // namespace laplace_likelihood

}  // namespace math
}  // namespace stan

#endif
