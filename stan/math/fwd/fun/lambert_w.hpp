#ifndef STAN_MATH_FWD_FUN_LAMBERT_W_HPP
#define STAN_MATH_FWD_FUN_LAMBERT_W_HPP

#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/meta.hpp>
#include <stan/math/prim/fun/lambert_w.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T>
inline fvar<T> lambert_w0(const fvar<T>& x) {
  const auto cached_result = lambert_w0(x.val_);
  return fvar<T>(cached_result, (x.d_ / (x.val_ + exp(cached_result))));
}

template <typename T>
inline fvar<T> lambert_wm1(const fvar<T>& x) {
  const auto cached_result = lambert_wm1(x.val_);
  return fvar<T>(cached_result, (x.d_ / (x.val_ + exp(cached_result))));
}

}  // namespace math
}  // namespace stan
#endif
