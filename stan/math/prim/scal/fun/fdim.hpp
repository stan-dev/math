#ifndef STAN_MATH_PRIM_SCAL_FUN_FDIM_HPP
#define STAN_MATH_PRIM_SCAL_FUN_FDIM_HPP

#include <stan/math/prim/meta.hpp>
#include <limits>

namespace stan {
namespace math {

/**
 * Return the positive difference of the specified values (C++11).
 *
 * The function is defined by
 *
 * <code>fdim(x, y) = (x > y) ? (x - y) : 0</code>.
 *
 * @param x First value.
 * @param y Second value.
 * @return max(x- y, 0)
 */
template <typename T1, typename T2, typename = require_all_arithmetic_t<T1, T2>>
inline double fdim(T1 x, T2 y) {
  using std::fdim;
  return fdim(x, y);
}

}  // namespace math
}  // namespace stan
#endif
