#ifndef STAN_MATH_PRIM_SCAL_FUN_FMIN_HPP
#define STAN_MATH_PRIM_SCAL_FUN_FMIN_HPP

#include <stan/math/prim/meta.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the lesser of the two specified arguments.  If one is
 * greater than the other, return not-a-number.
 *
 * @param x First argument.
 * @param y Second argument.
 * @return Minimum of x or y and if one is NaN return the other
 */
template <typename T1, typename T2, typename = require_all_arithmetic_t<T1, T2>>
inline double fmin(T1 x, T2 y) {
  using std::fmin;
  return fmin(x, y);
}

}  // namespace math
}  // namespace stan
#endif
