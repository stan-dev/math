#ifndef STAN_MATH_PRIM_FUN_ABS_HPP
#define STAN_MATH_PRIM_FUN_ABS_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/hypot.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return floating-point absolute value.
 *
 * Delegates to <code>fabs(double)</code> rather than
 * <code>std::abs(int)</code>.
 *
 * @param x scalar
 * @return absolute value of scalar
 */
inline double abs(double x) { return std::fabs(x); }

namespace internal {
/**
 * Return the absolute value of the complex argument.
 *
 * @tparam V value type of argument
 * @param[in] z argument
 * @return absolute value of the argument
 */
template <typename V>
inline V complex_abs(const std::complex<V>& z) {
  return hypot(z.real(), z.imag());
}
}  // namespace internal

}  // namespace math
}  // namespace stan

#endif
