#ifndef STAN_MATH_FWD_FUN_MULTIPLY_LOG_HPP
#define STAN_MATH_FWD_FUN_MULTIPLY_LOG_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/log.hpp>
#include <stan/math/fwd/fun/value_of_rec.hpp>
#include <stan/math/prim/fun/multiply_log.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T>
inline fvar<T> multiply_log(const fvar<T>& x1, const fvar<T>& x2) {
  if (value_of_rec(x1) == 0.0 && value_of_rec(x2) == 0.0) {
    return fvar<T>(0.0);
  }
  return fvar<T>(multiply_log(x1.val_, x2.val_),
                 x1.d_ * log(x2.val_) + x1.val_ * x2.d_ / x2.val_);
}

template <typename T>
inline fvar<T> multiply_log(double x1, const fvar<T>& x2) {
  if (x1 == 0.0 && value_of_rec(x2) == 0.0) {
    return fvar<T>(0.0);
  }
  return fvar<T>(multiply_log(x1, x2.val_), x1 * x2.d_ / x2.val_);
}

template <typename T>
inline fvar<T> multiply_log(const fvar<T>& x1, double x2) {
  if (value_of_rec(x1) == 0.0 && x2 == 0.0) {
    return fvar<T>(0.0);
  }
  return fvar<T>(multiply_log(x1.val_, x2), x1.d_ * log(x2));
}

}  // namespace math
}  // namespace stan
#endif
