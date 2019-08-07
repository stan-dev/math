#ifndef STAN_MATH_PRIM_SCAL_FUN_LAMBERTW_HPP
#define STAN_MATH_PRIM_SCAL_FUN_LAMBERTW_HPP

#include <stan/math/prim/meta.hpp>
#include <boost/math/special_functions/lambert_w.hpp>
#include <cmath>
#include <limits>

namespace stan {
namespace math {

template <typename T>
T lambert_w0(T x) {
  return boost::math::lambert_w0(x);
}

template <typename T>
T lambert_wm1(T x) {
  return boost::math::lambert_wm1(x);
}

}  // namespace math
}  // namespace stan
#endif
