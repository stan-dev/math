#ifndef STAN_MATH_FWD_FUN_ASIN_HPP
#define STAN_MATH_FWD_FUN_ASIN_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/fun/asin.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {

template <typename T>
inline fvar<T> asin(const fvar<T>& x) {
  using std::asin;
  using std::sqrt;
  return fvar<T>(asin(x.val_), x.d_ / sqrt(1 - square(x.val_)));
}

/**
 * Return the arc sine of the complex argument.
 *
 * @tparam T autodiff value type
 * @param[in] z argument
 * @return arc sine of the argument
 */
template <typename T>
inline std::complex<fvar<T>> asin(const std::complex<fvar<T>>& z) {
  return internal::complex_asin(z);
}

}  // namespace math
}  // namespace stan
#endif
