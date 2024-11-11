#ifndef STAN_MATH_FWD_FUN_LMULTIPLY_HPP
#define STAN_MATH_FWD_FUN_LMULTIPLY_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/log.hpp>
#include <stan/math/prim/fun/lmultiply.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T>
inline fvar<T> lmultiply(const fvar<T>& x1, const fvar<T>& x2) {
  return fvar<T>(lmultiply(x1.val_, x2.val_),
                 x1.d_ * log(x2.val_) + x1.val_ * x2.d_ / x2.val_);
}

template <typename T>
inline fvar<T> lmultiply(double x1, const fvar<T>& x2) {
  return fvar<T>(lmultiply(x1, x2.val_), x1 * x2.d_ / x2.val_);
}

template <typename T>
inline fvar<T> lmultiply(const fvar<T>& x1, double x2) {
  return fvar<T>(lmultiply(x1.val_, x2), x1.d_ * log(x2));
}
}  // namespace math
}  // namespace stan
#endif
