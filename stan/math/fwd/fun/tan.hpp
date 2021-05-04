#ifndef STAN_MATH_FWD_FUN_TAN_HPP
#define STAN_MATH_FWD_FUN_TAN_HPP

#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/meta.hpp>
#include <stan/math/prim/fun/tan.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {

template <typename T>
inline fvar<T> tan(const fvar<T>& x) {
  using std::cos;
  using std::tan;
  return fvar<T>(tan(x.val_), x.d_ / (cos(x.val_) * cos(x.val_)));
}

/**
 * Return the tangent of the complex argument.
 *
 * @param[in] z argument
 * @return tangent of the argument
 */
template <typename T>
inline std::complex<fvar<T>> tan(const std::complex<fvar<T>>& z) {
  return stan::math::internal::complex_tan(z);
}

}  // namespace math
}  // namespace stan
#endif
