#ifndef STAN_MATH_PRIM_FUN_STEP_HPP
#define STAN_MATH_PRIM_FUN_STEP_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>

namespace stan {
namespace math {

/**
 * The step, or Heaviside, function.
 *
 * The function is defined by
 *
 * <code>step(y) = (y < 0.0) ? 0 : 1</code>.
 *
   \f[
   \mbox{step}(x) =
   \begin{cases}
     0 & \mbox{if } x \le 0 \\
     1 & \mbox{if } x \ge 0  \\[6pt]
     0 & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @tparam T type of value
 * @param y value
 * @return zero if the value is less than zero, and one otherwise
 */
template <typename T, require_stan_scalar_t<T>* = nullptr>
inline T step(const T& y) {
  return y < 0.0 ? 0 : 1;
}

/**
 * Structure to wrap `step()` so it can be vectorized.
 */
struct step_fun {
  /**
   * Return zero if the value is less than zero and one otherwise
   *
   * @tparam T type of argument
   * @param y argument
   * @return step(y)
   */
  template <typename T>
  static inline T fun(const T& y) {
    return step(y);
  }
};

/**
 * Return the elementwise application of `step()` to
 * specified argument container.
 *
 * @tparam T type of container
 * @param x container
 * @return Elementwise step of members of container.
 */
template <typename T, require_container_t<T>* = nullptr>
inline auto step(const T& x) {
  return apply_scalar_unary<step_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan
#endif
