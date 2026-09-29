#ifndef STAN_MATH_REV_FUN_ERFCX_HPP
#define STAN_MATH_REV_FUN_ERFCX_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/erfcx.hpp>
#include <stan/math/prim/functor/apply_scalar_binary.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * The scaled complementary error function for variables.
 *
 * The derivative is
 *
 * \f$\frac{d}{dx}\mbox{erfcx}(x) = 2x\,\mbox{erfcx}(x) -
 \frac{2}{\sqrt{\pi}}\f$
 *
 * which reuses the function value, so no extra `exp` or `erfc` evaluation is
 * needed. That difference cancels for `x >= 4`, so `internal::erfcx_derivative`
 * takes the derivative from the tail rational there instead. Without that the
 * error reaches 2.55e+11 ulp at `x = 1e6`.
 *
   \f[
   \mbox{erfcx}(x) =
   \begin{cases}
     \exp(x^2)\operatorname{erfc}(x) & \mbox{if } -\infty\leq x \leq \infty
     \\[6pt] \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{erfcx}(x)}{\partial x} =
   \begin{cases}
     2x\operatorname{erfcx}(x) - \frac{2}{\sqrt{\pi}} & \mbox{if }
     -\infty\leq x\leq \infty \\[6pt] \textrm{NaN} & \mbox{if }
     x = \textrm{NaN}
   \end{cases}
   \f]
 *
 * One overload covers `var` and `var_value<Matrix>`. For a matrix the
 * derivative stays a lazy `binaryExpr`, so the reverse pass allocates
 * nothing.
 *
 * @tparam T a `var` or a `var_value` of a matrix type
 * @param a The variable.
 * @return Scaled complementary error function applied to the variable,
 *   elementwise for a matrix.
 */
template <
    typename T, require_var_t<T>* = nullptr,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr>
inline auto erfcx(T&& a) {
  auto val = to_arena(erfcx(a.val()));
  return make_callback_var(val, [a, val](auto& vi) mutable {
    // apply_scalar_binary returns a lazy binaryExpr for Eigen inputs, so
    // this allocates nothing, and calls the functor directly for scalars.
    auto deriv = apply_scalar_binary(
        [](double x, double v) { return internal::erfcx_derivative(x, v); },
        a.val(), val);
    as_array_or_scalar(a.adj())
        += as_array_or_scalar(vi.adj()) * as_array_or_scalar(deriv);
  });
}

}  // namespace math
}  // namespace stan
#endif
