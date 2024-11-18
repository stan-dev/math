#ifndef STAN_MATH_PRIM_FUN_ARG_HPP
#define STAN_MATH_PRIM_FUN_ARG_HPP

#include <stan/math/prim/fun/atan2.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {
/**
 * Return the phase angle of the complex argument.
 *
 * @tparam V value type of argument
 * @param[in] z argument
 * @return phase angle of the argument
 */
template <typename V, require_arithmetic_t<V>* = nullptr>
inline V arg(const std::complex<V>& z) {
  return std::arg(z);
}

namespace internal {
/**
 * Return the phase angle of the complex argument.
 *
 * @tparam V value type of argument
 * @param[in] z argument
 * @return phase angle of the argument
 */
template <typename V>
inline V complex_arg(const std::complex<V>& z) {
  return atan2(z.imag(), z.real());
}
}  // namespace internal
}  // namespace math
}  // namespace stan

#endif
