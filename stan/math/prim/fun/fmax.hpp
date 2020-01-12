#ifndef STAN_MATH_PRIM_FUN_FMAX_HPP
#define STAN_MATH_PRIM_FUN_FMAX_HPP

#include <stan/math/prim/meta.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the greater of the two specified arguments.  If one is
 * not-a-number, return the other.
 *
 * @param x First argument.
 * @param y Second argument.
 * @return maximum of x or y and if one is NaN return the other
 */
template <typename T1, typename T2, typename = require_all_arithmetic_t<T1, T2>>
inline double fmax(T1 x, T2 y) {
  using std::fmax;
  return fmax(x, y);
}

}  // namespace math
}  // namespace stan
#endif
