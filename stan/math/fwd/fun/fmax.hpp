#ifndef STAN_MATH_FWD_FUN_FMAX_HPP
#define STAN_MATH_FWD_FUN_FMAX_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/is_nan.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/fmax.hpp>

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
  if (unlikely(is_nan(x1))) {
    if (unlikely(is_nan(x2))) {
      return fvar<T>(NOT_A_NUMBER, NOT_A_NUMBER);
    }
    return x2;
  }
  if (unlikely(is_nan(x2))) {
    return x1;
  }
  return x1 > x2 ? x1 : x2;
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
    if (unlikely(is_nan(x2))) {
      return fvar<T>(NOT_A_NUMBER, NOT_A_NUMBER);
    }
    return x2;
  }
  if (unlikely(is_nan(x2))) {
    return fvar<T>(x1, 0);
  }
  return x1 > x2 ? fvar<T>(x1, 0) : x2;
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
  if (unlikely(is_nan(x1))) {
    if (unlikely(is_nan(x2))) {
      return fvar<T>(NOT_A_NUMBER, NOT_A_NUMBER);
    }
    return fvar<T>(x2, 0);
  }
  if (unlikely(is_nan(x2))) {
    return x1;
  }
  return x1 >= x2 ? x1 : fvar<T>(x2, 0);
}

}  // namespace math
}  // namespace stan
#endif
