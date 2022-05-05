#ifndef STAN_MATH_REV_FUN_SQUARE_HPP
#define STAN_MATH_REV_FUN_SQUARE_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>

namespace stan {
namespace math {

/**
 * Return the square of the input variable.
 *
 * <p>Using <code>square(x)</code> is more efficient
 * than using <code>x * x</code>.
 *
   \f[
   \mbox{square}(x) =
   \begin{cases}
     x^2 & \mbox{if } -\infty\leq x \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{square}(x)}{\partial x} =
   \begin{cases}
     2x & \mbox{if } -\infty\leq x\leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @param x Variable to square.
 * @return Square of variable.
 */
inline var square(const var& x) {
  return make_callback_var(square(x.val()), [x](auto& vi) mutable {
    x.adj() += vi.adj() * 2.0 * x.val();
  });
}

/**
 * Return the elementwise square of x
 *
 * @tparam T type of x
 * @param x argument
 * @return elementwise square of x
 */
template <typename T, require_rev_matrix_t<T>* = nullptr>
inline auto square(const T& x) {
  auto x_arena = to_arena(x);
  return make_callback_rev_matrix<T>(
      (x_arena.val().array().square()).matrix(), [x_arena](auto&& vi) mutable {
        x_arena.adj()
            += (2.0 * x_arena.val().array() * vi.adj().array()).matrix();
      });
}

}  // namespace math
}  // namespace stan
#endif
