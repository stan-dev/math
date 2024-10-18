
#ifndef STAN_MATH_VECTORIZED_FUN_SIGN_HPP
#define STAN_MATH_VECTORIZED_FUN_SIGN_HPP
#include <stan/math/prim/fun/sign.hpp>
namespace stan {
namespace math {

/**
 * Return the elementwise application of `sign()` to
 * specified argument container.
 *
 * @tparam T type of container
 * @param x container
 * @return Elementwise sign of members of container.
 */
template <typename T, require_container_t<T>* = nullptr>
inline auto sign(const T& x) {
  return apply_scalar_unary<sign_fun, T>::apply(x);
}


} // namespace math
} // namespace stan
#endif 

