#ifndef STAN_MATH_FWD_FUN_FMIN_HPP
#define STAN_MATH_FWD_FUN_FMIN_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/fmin.hpp>
#include <stan/math/prim/fun/is_nan.hpp>

namespace stan {
namespace math {

template <typename T>
inline fvar<T> fmin(const fvar<T>& x1, const fvar<T>& x2) {
  if (unlikely(is_nan(x1))) {
    if (unlikely(is_nan(x2))) {
      return fvar<T>(NOT_A_NUMBER, NOT_A_NUMBER);
    }
    return x2;
  }
  if (unlikely(is_nan(x2))) {
    return x1;
  }
  return x1 < x2 ? x1 : x2;
}

template <typename T>
inline fvar<T> fmin(double x1, const fvar<T>& x2) {
  if (unlikely(is_nan(x1))) {
    if (unlikely(is_nan(x2))) {
      return fvar<T>(NOT_A_NUMBER, NOT_A_NUMBER);
    }
    return x2;
  }
  if (unlikely(is_nan(x2))) {
    return fvar<T>(x1, 0);
  }
  return x2 <= x1 ? x2 : fvar<T>(x1, 0);
}

template <typename T>
inline fvar<T> fmin(const fvar<T>& x1, double x2) {
  if (unlikely(is_nan(x1))) {
    if (unlikely(is_nan(x2))) {
      return fvar<T>(NOT_A_NUMBER, NOT_A_NUMBER);
    }
    return fvar<T>(x2, 0);
  }
  if (unlikely(is_nan(x2))) {
    return x1;
  }
  return x1 <= x2 ? x1 : fvar<T>(x2, 0);
}

}  // namespace math
}  // namespace stan
#endif
