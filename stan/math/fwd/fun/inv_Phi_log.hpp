#ifndef STAN_MATH_FWD_FUN_INV_PHI_LOG_HPP
#define STAN_MATH_FWD_FUN_INV_PHI_LOG_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/inv_Phi_log.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T>
inline fvar<T> inv_Phi_log(const fvar<T>& p) {
  // using std::exp;
  const T xv = inv_Phi_log(p.val_);
  return fvar<T>(xv, exp(p.val_) * p.d_ / exp(-0.5 * square(xv)) * SQRT_TWO_PI);
}
}  // namespace math
}  // namespace stan
#endif
