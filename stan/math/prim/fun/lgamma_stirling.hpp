#ifndef STAN_MATH_PRIM_FUN_LGAMMA_STIRLING_HPP
#define STAN_MATH_PRIM_FUN_LGAMMA_STIRLING_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the Stirling approximation to the lgamma function.
 *
   \f[
   \mbox{lgamma_stirling}(x) =
    \frac{1}{2} \log(2\pi) + (x-\frac{1}{2})*\log(x) - x
   \f]
 *
 * @tparam T type of value
 * @param x value
 * @return Stirling's approximation to lgamma(x).
 */
template <typename T>
return_type_t<T> lgamma_stirling(const T x) {
  return HALF_LOG_TWO_PI + (x - 0.5) * log(x) - x;
}

}  // namespace math
}  // namespace stan

#endif
