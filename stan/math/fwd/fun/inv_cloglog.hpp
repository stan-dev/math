#ifndef STAN_MATH_FWD_FUN_INV_CLOGLOG_HPP
#define STAN_MATH_FWD_FUN_INV_CLOGLOG_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/exp.hpp>
#include <stan/math/prim/fun/inv_cloglog.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T>
inline fvar<T> inv_cloglog(const fvar<T>& x) {
  return fvar<T>(inv_cloglog(x.val_), x.d_ * exp(x.val_ - exp(x.val_)));
}
}  // namespace math
}  // namespace stan
#endif
