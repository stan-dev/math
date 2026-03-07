#ifndef STAN_MATH_MIX_FUNCTOR_LAPLACE_LIKELIHOOD_HPP
#define STAN_MATH_MIX_FUNCTOR_LAPLACE_LIKELIHOOD_HPP

#include <stan/math/mix/functor/hessian_block_diag.hpp>
#include <stan/math/mix/functor/conditional_copy_and_promote.hpp>
#include <stan/math/prim/functor.hpp>
#include <stan/math/prim/fun.hpp>
#include <type_traits>

namespace stan {
namespace math {

/**
 * functions to compute the log density, first, second,
 * and third-order derivatives for a likelihoood specified by the user.
 */
namespace laplace_likelihood {

namespace internal {

/**
 * Type trait to detect if a likelihood functor `F` provides a custom
 * `diff` method that computes the gradient and negative Hessian
 * analytically, avoiding the cost of embedded reverse-mode autodiff.
 *
 * A functor with a custom `diff` method should provide:
 *   auto diff(theta, hessian_block_size, args...) const
 * returning std::pair<gradient, sparse_hessian>.
 */
template <typename F, typename = void>
struct has_custom_diff : std::false_type {};

template <typename F>
struct has_custom_diff<F, std::void_t<decltype(std::declval<const F&>().diff(
                              std::declval<const Eigen::VectorXd&>(), 1))>>
    : std::true_type {};

template <typename F>
inline constexpr bool has_custom_diff_v = has_custom_diff<F>::value;

/**
 * Type trait to detect if a likelihood functor `F` provides a custom
 * `third_diff` method for the third derivative w.r.t. theta.
 */
template <typename F, typename = void>
struct has_custom_third_diff : std::false_type {};

template <typename F>
struct has_custom_third_diff<
    F, std::void_t<decltype(std::declval<const F&>().third_diff(
           std::declval<const Eigen::VectorXd&>()))>> : std::true_type {};

template <typename F>
inline constexpr bool has_custom_third_diff_v = has_custom_third_diff<F>::value;
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
 * Computes theta gradient `f` wrt `theta` and `args...`
 * @note If `Args` contains \ref var types then their adjoints will be
 * calculated as a side effect.
 * @tparam F A functor with `opertor()(Args&&...)` returning a scalar
 * @tparam Theta A class assignable to an Eigen vector type
 * @tparam Stream Type of stream for messages.
 * @tparam Args Type of variadic arguments.
 * @param f Log likelihood function.
 * @param theta Latent Gaussian model.
 * @param msgs Stream for messages.
 * @param args Variadic arguments for the likelihood function.
 */
template <typename F, typename Theta, typename Stream, typename... Args,
          require_eigen_vector_vt<std::is_arithmetic, Theta>* = nullptr>
inline auto theta_grad(F&& f, Theta&& theta, Stream* msgs, Args&&... args) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  nested_rev_autodiff nested;
  arena_t<Matrix<var, Dynamic, 1>> theta_var = theta;
  var f_var = f(theta_var, args..., msgs);
  grad(f_var.vi_);
  return theta_var.adj().eval();
}

/**
 * Computes likelihood argument gradient of `f`
 * @note If `Args` contains \ref var types then their adjoints will be
 * calculated as a side effect.
 * @tparam F A functor with `opertor()(Args&&...)` returning a scalar
 * @tparam Theta A class assignable to an Eigen vector type
 * @tparam Stream Type of stream for messages.
 * @tparam Args Type of variadic arguments.
 * @param f Log likelihood function.
 * @param theta Latent Gaussian model.
 * @param msgs Stream for messages.
 * @param args Variadic arguments for the likelihood function.
 */
template <typename F, typename Theta, typename Stream, typename... Args,
          require_eigen_vector_vt<std::is_arithmetic, Theta>* = nullptr>
inline void ll_arg_grad(F&& f, Theta&& theta, Stream* msgs, Args&&... args) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  nested_rev_autodiff nested;
  var f_var = f(theta, args..., msgs);
  grad(f_var.vi_);
}

/**
 * Computes negative diagonal Hessian of `f` wrt`theta` and `args...`
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
inline auto diagonal_hessian(F&& f, Theta&& theta, Stream* msgs,
                             Args&&... args) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  const Eigen::Index theta_size = theta.size();
  auto v = Eigen::VectorXd::Ones(theta_size);
  Eigen::VectorXd hessian_v = Eigen::VectorXd::Zero(theta_size);
  hessian_times_vector(f, hessian_v, std::forward<Theta>(theta), std::move(v),
                       value_of(args)..., msgs);
  return (-hessian_v).eval();
}

/**
 * Computes negative block diagonal Hessian of `f` wrt`theta` and `args...`
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
inline auto block_hessian(F&& f, Theta&& theta,
                          const Eigen::Index hessian_block_size, Stream* msgs,
                          Args&&... args) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  const Eigen::Index theta_size = theta.size();
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
    return (-hessian_theta).eval();
  } else {
    return (-hessian_block_diag(f, std::forward<Theta>(theta),
                                hessian_block_size, value_of(args)..., msgs))
        .eval();
  }
}

/**
 * Computes theta gradient and negative block diagonal Hessian of `f` wrt
 * `theta` and `args...`
 * @note If `Args` contains \ref var types then their adjoints will be
 * calculated as a side effect.
 * @note If `F` provides a custom `diff` method, it will be used instead
 * of the generic autodiff path for better performance.
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
  using F_t = std::decay_t<F>;
  if constexpr (has_custom_diff_v<F_t>) {
    // Use the functor's specialized analytic derivatives
    return f.diff(std::forward<Theta>(theta), hessian_block_size,
                  std::forward<Args>(args)...);
  } else {
    // Fall back to generic autodiff
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
      hessian_times_vector(f, hessian_v, std::forward<Theta>(theta),
                           std::move(v), value_of(args)..., msgs);
      Eigen::SparseMatrix<double> hessian_theta(theta_size, theta_size);
      hessian_theta.reserve(Eigen::VectorXi::Constant(theta_size, 1));
      for (Eigen::Index i = 0; i < theta_size; i++) {
        hessian_theta.insert(i, i) = hessian_v(i);
      }
      return std::make_pair(std::move(theta_gradient), (-hessian_theta).eval());
    } else {
      return std::make_pair(
          std::move(theta_gradient),
          (-hessian_block_diag(f, std::forward<Theta>(theta),
                               hessian_block_size, value_of(args)..., msgs))
              .eval());
    }
  }
}

/**
 * Compute third order derivative of `f` wrt `theta` and `args...`
 * @note If `Args` contains \ref var types then their adjoints will be
 * calculated as a side effect.
 * @note If `F` provides a custom `third_diff` method, it will be used
 * instead of the generic `fvar<fvar<var>>` autodiff path.
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
  using F_t = std::decay_t<F>;
  if constexpr (has_custom_third_diff_v<F_t>) {
    // Use the functor's specialized analytic third derivative
    return f.third_diff(std::forward<Theta>(theta),
                        std::forward<Args>(args)...);
  } else {
    // Fall back to generic fvar<fvar<var>> autodiff
    nested_rev_autodiff nested;
    const Eigen::Index theta_size = theta.size();
    arena_t<Eigen::Matrix<var, Eigen::Dynamic, 1>> theta_var
        = std::forward<Theta>(theta);
    arena_t<Eigen::Matrix<fvar<fvar<var>>, Eigen::Dynamic, 1>> theta_ffvar(
        theta_size);
    for (Eigen::Index i = 0; i < theta_size; ++i) {
      theta_ffvar(i) = fvar<fvar<var>>(fvar<var>(theta_var(i), 1.0), 1.0);
    }
    fvar<fvar<var>> ftheta_ffvar = f(theta_ffvar, args..., msgs);
    grad(ftheta_ffvar.d_.d_.vi_);
    return theta_var.adj().eval();
  }
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
  arena_t<Matrix<var, Dynamic, 1>> theta_var = std::forward<Theta>(theta);
  int n_blocks = theta_size / hessian_block_size;
  arena_t<VectorXd> v(theta_size);
  arena_t<VectorXd> w(theta_size);
  Matrix<fvar<fvar<var>>, Dynamic, 1> theta_ffvar(theta_size);
  auto shallow_copy_args
      = stan::math::internal::shallow_copy_vargs<fvar<fvar<var>>>(
          std::forward_as_tuple(args...));
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
  arena_t<Matrix<var, Dynamic, 1>> theta_var = std::forward<Theta>(theta);
  Matrix<fvar<var>, Dynamic, 1> theta_fvar(theta_size);
  for (Eigen::Index i = 0; i < theta_size; i++) {
    theta_fvar(i) = fvar<var>(theta_var(i), v(i));
  }
  auto shallow_copy_args = stan::math::internal::shallow_copy_vargs<fvar<var>>(
      std::forward_as_tuple(args...));
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
inline auto theta_grad(F&& f, Theta&& theta, TupleArgs&& ll_tup,
                       Stream* msgs = nullptr) {
  return apply(
      [](auto&& f, auto&& theta, auto&& msgs, auto&&... args) {
        return internal::theta_grad(std::forward<decltype(f)>(f),
                                    std::forward<decltype(theta)>(theta), msgs,
                                    std::forward<decltype(args)>(args)...);
      },
      std::forward<TupleArgs>(ll_tup), std::forward<F>(f),
      std::forward<Theta>(theta), msgs);
}

template <typename F, typename Theta, typename TupleArgs, typename Stream,
          require_eigen_vector_t<Theta>* = nullptr,
          require_tuple_t<TupleArgs>* = nullptr>
inline auto ll_arg_grad(F&& f, Theta&& theta, TupleArgs&& ll_tup,
                        Stream* msgs = nullptr) {
  return apply(
      [](auto&& f, auto&& theta, auto&& msgs, auto&&... args) {
        return internal::ll_arg_grad(std::forward<decltype(f)>(f),
                                     std::forward<decltype(theta)>(theta), msgs,
                                     std::forward<decltype(args)>(args)...);
      },
      std::forward<TupleArgs>(ll_tup), std::forward<F>(f),
      std::forward<Theta>(theta), msgs);
}

template <typename F, typename Theta, typename TupleArgs, typename Stream,
          require_eigen_vector_t<Theta>* = nullptr,
          require_tuple_t<TupleArgs>* = nullptr>
inline auto diagonal_hessian(F&& f, Theta&& theta, TupleArgs&& ll_tuple,
                             Stream* msgs) {
  return apply(
      [](auto&& f, auto&& theta, auto* msgs, auto&&... args) {
        return internal::diagonal_hessian(
            std::forward<decltype(f)>(f), std::forward<decltype(theta)>(theta),
            msgs, std::forward<decltype(args)>(args)...);
      },
      std::forward<TupleArgs>(ll_tuple), std::forward<F>(f),
      std::forward<Theta>(theta), msgs);
}

template <typename F, typename Theta, typename TupleArgs, typename Stream,
          require_eigen_vector_t<Theta>* = nullptr,
          require_tuple_t<TupleArgs>* = nullptr>
inline auto block_hessian(F&& f, Theta&& theta,
                          const Eigen::Index hessian_block_size,
                          TupleArgs&& ll_tuple, Stream* msgs) {
  return apply(
      [](auto&& f, auto&& theta, auto hessian_block_size, auto* msgs,
         auto&&... args) {
        return internal::block_hessian(
            std::forward<decltype(f)>(f), std::forward<decltype(theta)>(theta),
            hessian_block_size, msgs, std::forward<decltype(args)>(args)...);
      },
      std::forward<TupleArgs>(ll_tuple), std::forward<F>(f),
      std::forward<Theta>(theta), hessian_block_size, msgs);
}

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
