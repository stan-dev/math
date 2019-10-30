#ifndef STAN_MATH_PRIM_SCAL_FUN_HYPOT_HPP
#define STAN_MATH_PRIM_SCAL_FUN_HYPOT_HPP

#include <stan/math/prim/meta.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the length of the hypoteneuse of a right triangle with
 * opposite and adjacent side lengths given by the specified
 * arguments (C++11).  In symbols, if the arguments are
 * <code>x</code> and <code>y</code>, the result is <code>sqrt(x *
 * x + y * y)</code>.
 *
 * @param x First argument.
 * @param y Second argument.
 * @return Length of hypoteneuse of right triangle with opposite
 * and adjacent side lengths x and y.
 */
template <typename T1, typename T2, typename = require_all_arithmetic_t<T1, T2>>
inline double hypot(T1 x, T2 y) {
  using std::hypot;
  return hypot(x, y);
}

}  // namespace math
}  // namespace stan
#endif
