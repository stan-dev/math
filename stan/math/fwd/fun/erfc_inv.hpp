#ifndef STAN_MATH_FWD_FUN_ERFC_INV_HPP
#define STAN_MATH_FWD_FUN_ERFC_INV_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/erfc_inv.hpp>
#include <stan/math/fwd/fun/square.hpp>
#include <stan/math/fwd/fun/exp.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T>
inline fvar<T> erfc_inv(const fvar<T>& x) {
  T precomp_erfc_inv = erfc_inv(x.val());
  return fvar<T>(precomp_erfc_inv,
                 -x.d_ * exp(LOG_SQRT_PI - LOG_TWO + square(precomp_erfc_inv)));
}

}  // namespace math
}  // namespace stan
#endif
