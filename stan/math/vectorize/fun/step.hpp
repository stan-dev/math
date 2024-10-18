
#ifndef STAN_MATH_VECTORIZED_FUN_STEP_HPP
#define STAN_MATH_VECTORIZED_FUN_STEP_HPP
#include <stan/math/prim/fun/step.hpp>
namespace stan {
namespace math {

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


} // namespace math
} // namespace stan
#endif 

