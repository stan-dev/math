#ifndef STAN_MATH_FWD_FUN_ABS_HPP
#define STAN_MATH_FWD_FUN_ABS_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/abs.hpp>
#include <stan/math/prim/fun/value_of.hpp>

namespace stan {
namespace math {

template <typename T>
inline fvar<T> abs(const fvar<T>& x) {
  if (x.val_ > 0.0) {
    return x;
  } else if (x.val_ < 0.0) {
    return fvar<T>(-x.val_, -x.d_);
  } else if (x.val_ == 0.0) {
    return fvar<T>(0, 0);
  } else {
    return fvar<T>(abs(x.val_), NOT_A_NUMBER);
  }
}

}  // namespace math
}  // namespace stan
#endif
