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
inline auto inv(const var a) {
  auto denom = to_arena(square(a.val()));
  return make_callback_var(inv(a.val()), [a, denom](auto& vi) mutable {
    a.adj() -= vi.adj() / denom;
  });
}

template <typename T, require_rev_matrix_t<T>* = nullptr>
inline auto inv(const T& x) {
  auto x_arena = to_arena(x);
  auto denom = to_arena(square(x_arena.val()).array());
  return make_callback_rev_matrix<T>(
      inv(x_arena.val()), [x_arena, denom](auto&& vi) mutable {
        x_arena.adj().array() -= vi.adj().array() / denom;
      });
}
}  // namespace math
}  // namespace stan
#endif
