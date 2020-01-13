#ifndef STAN_MATH_PRIM_FUN_LGAMMA_STIRLING_HPP
#define STAN_MATH_PRIM_FUN_LGAMMA_STIRLING_HPP

#include <cmath>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/lgamma.hpp>

namespace stan {
namespace math {

/**
 * Return the Stirling approximation to the gamma function.
 *

   \f[
   \mbox{lgamma_stirling}(x) =
    \frac{1}{2} \log(2\pi) + (x-\frac{1}{2})*\log(x) - x
   \f]

 *
 * @param x value
 * @return Stirling's approximation to lgamma(x).
 * @tparam T Type of  value.
 */
template <typename T>
T lgamma_stirling(const T x) {
  using std::log;
  return 0.5 * log(2 * stan::math::pi()) + (x - 0.5) * log(x) - x;
}

}  // namespace math
}  // namespace stan

#endif
