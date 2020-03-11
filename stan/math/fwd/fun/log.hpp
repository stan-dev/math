#ifndef STAN_MATH_FWD_FUN_LOG_HPP
#define STAN_MATH_FWD_FUN_LOG_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {

template <typename T>
inline fvar<T> log(const fvar<T>& x) {
  using std::log;
  if (x.val_ < 0.0) {
    return fvar<T>(NOT_A_NUMBER, NOT_A_NUMBER);
  } else {
    return fvar<T>(log(x.val_), x.d_ / x.val_);
  }
}

/**
 * Return the natural logarithm (base e) of the specified complex argument.
 *
 * @tparam T autodiff value type
 * @param z complex argument
 * @return natural logarithm of argument
 */
template <typename T>
inline std::complex<fvar<T>> log(const std::complex<fvar<T>>& z) {
  return internal::complex_log(z);
}

}  // namespace math
}  // namespace stan
#endif
