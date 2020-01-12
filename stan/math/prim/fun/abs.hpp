#ifndef STAN_MATH_PRIM_FUN_ABS_HPP
#define STAN_MATH_PRIM_FUN_ABS_HPP

#include <stan/math/prim/meta.hpp>
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

}  // namespace math
}  // namespace stan

#endif
