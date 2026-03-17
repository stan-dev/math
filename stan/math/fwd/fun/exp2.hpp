#ifndef STAN_MATH_FWD_FUN_EXP2_HPP
#define STAN_MATH_FWD_FUN_EXP2_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/exp2.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T>
inline fvar<T> exp2(const fvar<T>& x) {
  return fvar<T>(exp2(x.val_), x.d_ * exp2(x.val_) * LOG_TWO);
}

}  // namespace math
}  // namespace stan
#endif
