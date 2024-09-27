#ifndef STAN_MATH_FWD_FUN_SIN_HPP
#define STAN_MATH_FWD_FUN_SIN_HPP

#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/meta.hpp>
#include <stan/math/prim/fun/sin.hpp>
#include <cmath>
#include <stan/math/prim/core/complex_base.hpp>

namespace stan {
namespace math {

template <typename T>
inline fvar<T> sin(const fvar<T>& x) {
  using std::cos;
  using std::sin;
  return fvar<T>(sin(x.val_), x.d_ * cos(x.val_));
}

/**
 * Return the sine of the complex argument.
 *
 * @tparam T autodiff value type
 * @param[in] z argument
 * @return sine of the argument
 */
template <typename T>
inline stan::math::complex<fvar<T>> sin(const stan::math::complex<fvar<T>>& z) {
  return internal::complex_sin(z);
}

}  // namespace math
}  // namespace stan
#endif
