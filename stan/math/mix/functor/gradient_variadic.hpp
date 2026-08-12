#ifndef STAN_MATH_MIX_FUNCTOR_GRADIENT_VARIADIC_HPP
#define STAN_MATH_MIX_FUNCTOR_GRADIENT_VARIADIC_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/apply.hpp>
#include <stan/math/prim/functor/conditional_copy_and_promote.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/to_arena.hpp>
#include <stan/math/fwd/core.hpp>
#include <ostream>
#include <tuple>

namespace stan {
namespace math {

/**
 * Return the gradient of a scalar functor with respect to its first argument.
 *
 * The functor must be callable as `f(theta, msgs, args...)`. This overload is
 * used when none of the inputs contain reverse-mode variables.
 *
 * @tparam F type of functor
 * @tparam ThetaVec Eigen column vector type
 * @tparam Args types of additional arguments
 * @param f functor returning a scalar
 * @param theta argument with respect to which the gradient is taken
 * @param msgs output stream for functor messages
 * @param args additional arguments passed to the functor
 * @return gradient of `f` with respect to `theta`
 */
template <typename F, typename ThetaVec, typename... Args,
          require_eigen_col_vector_t<ThetaVec>* = nullptr,
          require_t<is_vt_not_complex<ThetaVec>>* = nullptr,
          require_not_t<is_any_var_scalar<ThetaVec, Args...>>* = nullptr>
inline Eigen::VectorXd gradient(const F& f, const ThetaVec& theta,
                                std::ostream* msgs, const Args&... args) {
  const Eigen::Index n = theta.size();
  Eigen::VectorXd g(n);
  if (n == 0) {
    return g;
  }

  nested_rev_autodiff nested;
  Eigen::Matrix<var, Eigen::Dynamic, 1> theta_var = theta;
  var fx = f(theta_var, msgs, args...);
  grad(fx.vi_);
  g = theta_var.adj();
  return g;
}

/**
 * Return a differentiable gradient of a scalar functor with respect to its
 * first argument.
 *
 * For output adjoints `w`, the reverse pass propagates
 * `H(theta) * w` to `theta` and the corresponding mixed partials to
 * `args...`. The functor must be callable with `fvar<var>` inputs.
 *
 * @tparam F type of functor
 * @tparam ThetaVec Eigen column vector type
 * @tparam Args types of additional arguments
 * @param f functor returning a scalar
 * @param theta argument with respect to which the gradient is taken
 * @param msgs output stream for functor messages
 * @param args additional arguments passed to the functor
 * @return gradient of `f` with respect to `theta`
 */
template <typename F, typename ThetaVec, typename... Args,
          require_eigen_col_vector_t<ThetaVec>* = nullptr,
          require_t<is_vt_not_complex<ThetaVec>>* = nullptr,
          require_t<is_any_var_scalar<ThetaVec, Args...>>* = nullptr>
inline Eigen::Matrix<var, Eigen::Dynamic, 1> gradient(
    const F& f, const ThetaVec& theta, std::ostream* msgs,
    const Args&... args) {
  using Eigen::Index;
  const Index n = theta.size();
  if (n == 0) {
    return Eigen::Matrix<var, Eigen::Dynamic, 1>(0);
  }

  auto theta_arena = to_arena(theta);
  auto* args_arena = make_chainable_ptr(std::make_tuple(args...));
  Eigen::VectorXd g(n);
  {
    nested_rev_autodiff nested;
    Eigen::Matrix<var, Eigen::Dynamic, 1> theta_var
        = value_of(theta_arena);
    var fx = math::apply(
        [&](const auto&... inner_args) {
          return f(theta_var, msgs, value_of(inner_args)...);
        },
        *args_arena);
    grad(fx.vi_);
    g = theta_var.adj();
  }

  arena_t<Eigen::Matrix<var, Eigen::Dynamic, 1>> res = g;
  reverse_pass_callback([f, msgs, theta_arena, args_arena, res]() mutable {
    const Index n = theta_arena.size();
    nested_rev_autodiff nested;
    Eigen::Matrix<fvar<var>, Eigen::Dynamic, 1> theta_fvar(n);
    for (Index i = 0; i < n; ++i) {
      theta_fvar.coeffRef(i)
          = fvar<var>(theta_arena.coeff(i), res.adj().coeff(i));
    }
    auto args_fvar
        = internal::shallow_copy_vargs<fvar<var>>(*args_arena);
    fvar<var> fx = math::apply(
        [&](const auto&... inner_args) {
          return f(theta_fvar, msgs, inner_args...);
        },
        args_fvar);
    grad(fx.d_.vi_);
  });
  return res;
}

}  // namespace math
}  // namespace stan
#endif
