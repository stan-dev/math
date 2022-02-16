#ifndef STAN_MATH_FWD_FUN_INV_ERFC_HPP
#define STAN_MATH_FWD_FUN_INV_ERFC_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/inv_erfc.hpp>
#include <stan/math/fwd/fun/square.hpp>
#include <stan/math/fwd/fun/exp.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T>
inline fvar<T> inv_erfc(const fvar<T>& x) {
  T precomp_inv_erfc = inv_erfc(x.val());
  return fvar<T>(precomp_inv_erfc,
                 -x.d_ * exp(LOG_SQRT_PI - LOG_TWO + square(precomp_inv_erfc)));
}

}  // namespace math
}  // namespace stan
#endif
