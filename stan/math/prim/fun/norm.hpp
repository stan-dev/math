#ifndef STAN_MATH_PRIM_FUN_NORM_HPP
#define STAN_MATH_PRIM_FUN_NORM_HPP

#include <stan/math/prim/fun/isinf.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {
namespace internal {
/**
 * Return the squared magnitude of the complex argument.
 *
 * @tparam V value type of argument
 * @param[in] z argument
 * @return squared magnitude of the argument
 */
template <typename V>
inline V complex_norm(const std::complex<V>& z) {
  return square(z.real()) + square(z.imag());
}
}  // namespace internal
}  // namespace math
}  // namespace stan

#endif
