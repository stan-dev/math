#ifndef STAN_MATH_FWD_FUN_LOG10_HPP
#define STAN_MATH_FWD_FUN_LOG10_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T>
inline fvar<T> log10(const fvar<T>& x) {
  using std::log;
  using std::log10;
  if (x.val_ < 0.0) {
    return fvar<T>(NOT_A_NUMBER, NOT_A_NUMBER);
  } else {
    return fvar<T>(log10(x.val_), x.d_ / (x.val_ * LOG_TEN));
  }
}

}  // namespace math
}  // namespace stan
#endif
