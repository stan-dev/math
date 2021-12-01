#ifndef STAN_MATH_REV_FUN_INV_HPP
#define STAN_MATH_REV_FUN_INV_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/inv.hpp>

namespace stan {
namespace math {

/**
 * @tparam T Arithmetic or a type inheriting from `EigenBase`.
   \f[
   \mbox{inv}(x) =
   \begin{cases}
     \frac{1}{x} & \mbox{if } -\infty\leq x \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{inv}(x)}{\partial x} =
   \begin{cases}
     -\frac{1}{x^2} & \mbox{if } -\infty\leq x\leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]
 *
 */
template <typename T, require_stan_scalar_or_eigen_t<T>* = nullptr>
inline auto inv(const var_value<T>& a) {
  auto denom = to_arena(as_array_or_scalar(square(a.val())));
  return make_callback_var(inv(a.val()), [a, denom](auto& vi) mutable {
    as_array_or_scalar(a.adj()) -= as_array_or_scalar(vi.adj()) / denom;
  });
}

}  // namespace math
}  // namespace stan
#endif
