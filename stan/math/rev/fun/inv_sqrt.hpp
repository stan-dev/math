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
inline auto inv_sqrt(const var a) {
  auto denom = to_arena(a.val() * sqrt(a.val()));
  return make_callback_var(inv_sqrt(a.val()), [a, denom](auto& vi) mutable {
    a.adj() -= 0.5 * vi.adj() / denom;
  });
}

template <typename T, require_rev_matrix_t<T>* = nullptr>
inline auto inv_sqrt(const T& x) {
  auto x_arena = to_arena(x);
  auto denom = to_arena(x_arena.val().array() * sqrt(x_arena.val()).array());
  return make_callback_rev_matrix<T>(inv_sqrt(x_arena.val()), [x_arena, denom](auto&& vi) mutable {
    x_arena.adj().array() -= 0.5 * vi.adj().array() / denom;
  });
}

}  // namespace math
}  // namespace stan
#endif
