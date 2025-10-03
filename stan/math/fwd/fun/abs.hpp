#ifndef STAN_MATH_FWD_FUN_ABS_HPP
#define STAN_MATH_FWD_FUN_ABS_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/value_of.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/abs.hpp>
#include <complex>

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

/**
 * Return the absolute value of the complex argument.
 *
 * @tparam T value type of argument
 * @param[in] z argument
 * @return absolute value of the argument
 */
template <typename T>
inline fvar<T> abs(const std::complex<fvar<T>>& z) {
  return internal::complex_abs(z);
}

}  // namespace math
}  // namespace stan
#endif
