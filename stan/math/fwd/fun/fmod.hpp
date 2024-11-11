#ifndef STAN_MATH_FWD_FUN_FMOD_HPP
#define STAN_MATH_FWD_FUN_FMOD_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/trunc.hpp>
#include <stan/math/fwd/fun/value_of.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/is_any_nan.hpp>
#include <stan/math/prim/fun/fmod.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T>
inline fvar<T> fmod(const fvar<T>& x1, const fvar<T>& x2) {
  return fvar<T>(fmod(x1.val_, x2.val_),
                 x1.d_ - x2.d_ * trunc(x1.val_ / x2.val_));
}

template <typename T>
inline fvar<T> fmod(const fvar<T>& x1, double x2) {
  if (unlikely(is_any_nan(value_of(x1.val_), x2))) {
    return fvar<T>(fmod(x1.val_, x2), NOT_A_NUMBER);
  } else {
    return fvar<T>(fmod(x1.val_, x2), x1.d_);
  }
}

template <typename T>
inline fvar<T> fmod(double x1, const fvar<T>& x2) {
  return fvar<T>(fmod(x1, x2.val_), -x2.d_ * trunc(x1 / x2.val_));
}

}  // namespace math
}  // namespace stan
#endif
