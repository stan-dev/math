#ifndef STAN_MATH_FWD_FUN_ERFC_HPP
#define STAN_MATH_FWD_FUN_ERFC_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/exp.hpp>
#include <stan/math/fwd/fun/square.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/erfc.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T>
inline fvar<T> erfc(const fvar<T>& x) {
  return fvar<T>(erfc(x.val_), -x.d_ * exp(-square(x.val_)) * TWO_OVER_SQRT_PI);
}

}  // namespace math
}  // namespace stan
#endif
