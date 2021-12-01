#ifndef STAN_MATH_REV_FUN_INV_SQRT_HPP
#define STAN_MATH_REV_FUN_INV_SQRT_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/inv_sqrt.hpp>
#include <stan/math/rev/fun/to_arena.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * @tparam T Arithmetic or a type inheriting from `EigenBase`.
 *
   \f[
   \mbox{inv\_sqrt}(x) =
   \begin{cases}
     \frac{1}{\sqrt{x}} & \mbox{if } -\infty\leq x \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{inv\_sqrt}(x)}{\partial x} =
   \begin{cases}
     -\frac{1}{2\sqrt{x^3}} & \mbox{if } -\infty\leq x\leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]
 *
 */
template <typename T, require_stan_scalar_or_eigen_t<T>* = nullptr>
inline auto inv_sqrt(const var_value<T>& a) {
  auto denom = to_arena(as_array_or_scalar(a.val())
                        * as_array_or_scalar(sqrt(a.val())));
  return make_callback_var(inv_sqrt(a.val()), [a, denom](auto& vi) mutable {
    as_array_or_scalar(a.adj()) -= 0.5 * as_array_or_scalar(vi.adj()) / denom;
  });
}

}  // namespace math
}  // namespace stan
#endif
