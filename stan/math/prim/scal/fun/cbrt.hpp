#ifndef STAN_MATH_PRIM_SCAL_FUN_CBRT_HPP
#define STAN_MATH_PRIM_SCAL_FUN_CBRT_HPP

#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/fun/is_inf.hpp>
#include <stan/math/prim/scal/fun/is_nan.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the cube root of the specified value
 *
 * @param[in] x Argument.
 * @return Cube root of the argument.
 * @throw std::domain_error If argument is negative.
 */
inline double cbrt(double x) {
  if (is_nan(x))
    return NOT_A_NUMBER;
  if (is_inf(x))
    return x < 0 ? NEGATIVE_INFTY : INFTY;
  return std::cbrt(x);
}

/**
 * Integer version of cbrt.
 *
 * @param[in] x Argument.
 * @return Cube root of the argument.
 * @throw std::domain_error If argument is less than 1.
 */
inline double cbrt(int x) { return cbrt(static_cast<double>(x)); }

}  // namespace math
}  // namespace stan
#endif
