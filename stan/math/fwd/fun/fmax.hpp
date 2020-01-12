#ifndef STAN_MATH_FWD_FUN_FMAX_HPP
#define STAN_MATH_FWD_FUN_FMAX_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/fmax.hpp>
#include <stan/math/prim/fun/is_nan.hpp>

namespace stan {
namespace math {

/**
 * Return the greater of the two specified arguments.  If one is
 * not-a-number, return the other.
 *
 * @tparam T inner type of the fvar
 * @param x1 First argument.
 * @param x2 Second argument.
 * @return maximum of arguments, and if one is NaN return the other
 */
template <typename T>
inline fvar<T> fmax(const fvar<T>& x1, const fvar<T>& x2) {
  if (unlikely(is_nan(x1.val_))) {
    if (is_nan(x2.val_)) {
      return fvar<T>(fmax(x1.val_, x2.val_), NOT_A_NUMBER);
    } else {
      return fvar<T>(x2.val_, x2.d_);
    }
  } else if (unlikely(is_nan(x2.val_))) {
    return fvar<T>(x1.val_, x1.d_);
  } else if (x1.val_ > x2.val_) {
    return fvar<T>(x1.val_, x1.d_);
  } else if (x1.val_ == x2.val_) {
    return fvar<T>(x1.val_, NOT_A_NUMBER);
  } else {
    return fvar<T>(x2.val_, x2.d_);
  }
}

/**
 * Return the greater of the two specified arguments.  If one is
 * not-a-number, return the other.
 *
 * @tparam T inner type of the fvar
 * @param x1 First argument.
 * @param x2 Second argument.
 * @return maximum of arguments, and if one is NaN return the other
 */
template <typename T>
inline fvar<T> fmax(double x1, const fvar<T>& x2) {
  if (unlikely(is_nan(x1))) {
    if (is_nan(x2.val_)) {
      return fvar<T>(fmax(x1, x2.val_), NOT_A_NUMBER);
    } else {
      return fvar<T>(x2.val_, x2.d_);
    }
  } else if (unlikely(is_nan(x2.val_))) {
    return fvar<T>(x1, 0.0);
  } else if (x1 > x2.val_) {
    return fvar<T>(x1, 0.0);
  } else if (x1 == x2.val_) {
    return fvar<T>(x2.val_, NOT_A_NUMBER);
  } else {
    return fvar<T>(x2.val_, x2.d_);
  }
}

/**
 * Return the greater of the two specified arguments.  If one is
 * not-a-number, return the other.
 *
 * @tparam T inner type of the fvar
 * @param x1 First argument.
 * @param x2 Second argument.
 * @return maximum of arguments, and if one is NaN return the other
 */
template <typename T>
inline fvar<T> fmax(const fvar<T>& x1, double x2) {
  if (unlikely(is_nan(x1.val_))) {
    if (is_nan(x2)) {
      return fvar<T>(fmax(x1.val_, x2), NOT_A_NUMBER);
    } else {
      return fvar<T>(x2, 0.0);
    }
  } else if (unlikely(is_nan(x2))) {
    return fvar<T>(x1.val_, x1.d_);
  } else if (x1.val_ > x2) {
    return fvar<T>(x1.val_, x1.d_);
  } else if (x1.val_ == x2) {
    return fvar<T>(x1.val_, NOT_A_NUMBER);
  } else {
    return fvar<T>(x2, 0.0);
  }
}

}  // namespace math
}  // namespace stan
#endif
