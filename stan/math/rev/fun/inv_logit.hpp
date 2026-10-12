#ifndef STAN_MATH_REV_FUN_INV_LOGIT_HPP
#define STAN_MATH_REV_FUN_INV_LOGIT_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/abs.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/inv_logit.hpp>
#include <stan/math/prim/functor/apply_scalar_binary.hpp>

namespace stan {
namespace math {

/**
 * The inverse logit function for variables (stan).
 *
 * See inv_logit() for the double-based version.
 *
 * The derivative of inverse logit is
 *
 * \f$\frac{d}{dx} \mbox{logit}^{-1}(x) = \frac{e}{(1 + e)^2}\f$
 * with \f$e = \exp(-|x|)\f$.
 *
 * @tparam T a `var` or a `var_value` of a matrix type
 * @param a Argument variable.
 * @return Inverse logit of argument, elementwise for a matrix.
 */
template <
    typename T, require_var_t<T>* = nullptr,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr>
inline auto inv_logit(T&& a) {
  auto e = to_arena(exp(-abs(a.val())));
  return make_callback_var(
      apply_scalar_binary(
          [](double x, double e) { return (x >= 0 ? 1.0 : e) / (1.0 + e); },
          a.val(), e),
      [a, e](auto& vi) mutable {
        const auto& e_a = as_array_or_scalar(e);
        as_array_or_scalar(a.adj())
            += as_array_or_scalar(vi.adj()) * e_a / ((1.0 + e_a) * (1.0 + e_a));
      });
}

/**
 * The inverse logit function for Eigen expressions with var value type.
 *
 * See inv_logit() for the double-based version.
 *
 * The derivative of inverse logit is
 *
 * \f$\frac{d}{dx} \mbox{logit}^{-1}(x) = \frac{e}{(1 + e)^2}\f$
 * with \f$e = \exp(-|x|)\f$.
 *
 * @tparam T type of Eigen expression
 * @param x Eigen expression
 * @return Inverse logit of each value in x.
 */
template <typename T, require_eigen_vt<is_var, T>* = nullptr>
inline auto inv_logit(T&& x) {
  auto x_arena = to_arena(std::forward<T>(x));
  const auto x_val = x_arena.val().array().eval();
  arena_t<plain_type_t<decltype(x_val)>> e = (-x_val.abs()).exp();
  arena_t<T> ret = ((x_val >= 0).select(1.0, e) / (1.0 + e)).matrix();
  reverse_pass_callback([x_arena, ret, e]() mutable {
    x_arena.adj().array() += ret.adj().array() * e / (1.0 + e).square();
  });
  return ret;
}

}  // namespace math
}  // namespace stan
#endif
