#ifndef STAN_MATH_PRIM_SCAL_FUN_EXP2_HPP
#define STAN_MATH_PRIM_SCAL_FUN_EXP2_HPP

#include <stan/math/prim/meta.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the exponent base 2 of the specified argument (C99,
 * C++11).
 *
 * The exponent base 2 function is defined by
 *
 * <code>exp2(y) = pow(2.0, y)</code>.
 *
 * @param y argument.
 * @return exponent base 2 of argument.
 */
template <typename T, typename = require_arithmetic_t<T>>
inline double exp2(T y) {
  using std::exp2;
  return exp2(y);
}

}  // namespace math
}  // namespace stan
#endif
