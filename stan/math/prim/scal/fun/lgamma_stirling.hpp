#ifndef STAN_MATH_PRIM_SCAL_FUN_LGAMMA_STIRLING_HPP
#define STAN_MATH_PRIM_SCAL_FUN_LGAMMA_STIRLING_HPP

#include <cmath>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/fun/lgamma.hpp>

namespace stan {
namespace math {

template <typename T>
T lgamma_stirling(const T x) {
  return 0.5 * log(2 * stan::math::pi()) + (x - 0.5) * log(x) - x;
}

}  // namespace math
}  // namespace stan

#endif
