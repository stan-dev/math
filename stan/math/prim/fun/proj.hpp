#ifndef STAN_MATH_PRIM_FUN_PROJ_HPP
#define STAN_MATH_PRIM_FUN_PROJ_HPP

#include <stan/math/prim/fun/is_inf.hpp>
#include <complex>
#include <limits>

namespace stan {
namespace math {
/**
 * Return the projection of the complex argument onto the Riemann
 * sphere.
 *
 * @tparam V value type of argument
 * @param[in] z argument
 * @return projection of the argument onto the Riemann sphere
 */
template <typename V>
inline std::complex<V> proj(const std::complex<V>& z) {
  return std::proj(z);
}

namespace internal {
/**
 * Return the projection of the complex argument onto the Riemann
 * sphere.
 *
 * @tparam V value type of argument
 * @param[in] z argument
 * @return projection of the argument onto the Riemann sphere
 */
template <typename V>
inline std::complex<V> complex_proj(const std::complex<V>& z) {
  if (is_inf(z.real()) || is_inf(z.imag())) {
    return {std::numeric_limits<V>::infinity(), z.imag() < 0 ? -0.0 : 0.0};
  }
  return z;
}
}  // namespace internal
}  // namespace math
}  // namespace stan

#endif
