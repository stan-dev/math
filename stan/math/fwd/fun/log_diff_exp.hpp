#ifndef STAN_MATH_FWD_FUN_LOG_DIFF_EXP_HPP
#define STAN_MATH_FWD_FUN_LOG_DIFF_EXP_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/expm1.hpp>
#include <stan/math/prim/fun/log_diff_exp.hpp>

namespace stan {
namespace math {

template <typename T>
inline fvar<T> log_diff_exp(const fvar<T>& x1, const fvar<T>& x2) {
  if (x1.val_ <= x2.val_) {
    if (x1.val_ < INFTY && x1.val_ == x2.val_) {
      return fvar<T>(NEGATIVE_INFTY, NOT_A_NUMBER);
    }
    return fvar<T>(NOT_A_NUMBER, NOT_A_NUMBER);
  }
  return fvar<T>(
      log_diff_exp(x1.val_, x2.val_),
      -(x1.d_ / expm1(x2.val_ - x1.val_) + x2.d_ / expm1(x1.val_ - x2.val_)));
}

template <typename T1, typename T2, require_arithmetic_t<T1>* = nullptr>
inline fvar<T2> log_diff_exp(const T1& x1, const fvar<T2>& x2) {
  if (x1 <= x2.val_) {
    if (x1 < INFTY && x1 == x2.val_) {
      return fvar<T2>(NEGATIVE_INFTY, x2.d_ * NEGATIVE_INFTY);
    }
    return fvar<T2>(NOT_A_NUMBER, NOT_A_NUMBER);
  }
  return fvar<T2>(log_diff_exp(x1, x2.val_), -x2.d_ / expm1(x1 - x2.val_));
}

template <typename T1, typename T2, require_arithmetic_t<T2>* = nullptr>
inline fvar<T1> log_diff_exp(const fvar<T1>& x1, const T2& x2) {
  if (x1.val_ <= x2) {
    if (x1.val_ < INFTY && x1.val_ == x2) {
      if (x2 == NEGATIVE_INFTY) {
        return fvar<T1>(NEGATIVE_INFTY, x1.d_);
      }
      return fvar<T1>(NEGATIVE_INFTY, x1.d_ * INFTY);
    }
    return fvar<T1>(NOT_A_NUMBER, NOT_A_NUMBER);
  }
  return fvar<T1>(log_diff_exp(x1.val_, x2), -x1.d_ / expm1(x2 - x1.val_));
}
}  // namespace math
}  // namespace stan
#endif
