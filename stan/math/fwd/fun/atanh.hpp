#ifndef STAN_MATH_FWD_FUN_ATANH_HPP
#define STAN_MATH_FWD_FUN_ATANH_HPP

#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/meta.hpp>
#include <stan/math/prim/fun/atanh.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <cmath>
#include <stan/math/prim/core/complex_base.hpp>

namespace stan {
namespace math {

/**
 * Return inverse hyperbolic tangent of specified value.
 *
 * @tparam T scalar type of forward-mode autodiff variable
 * argument.
 * @param x Argument.
 * @return Inverse hyperbolic tangent of argument.
 * @throw std::domain_error if x < -1 or x > 1.
 */
template <typename T>
inline fvar<T> atanh(const fvar<T>& x) {
  return fvar<T>(atanh(x.val_), x.d_ / (1 - square(x.val_)));
}

/**
 * Return the hyperbolic arc tangent of the complex argument.
 *
 * @tparam T autodiff value type
 * @param[in] z argument
 * @return hyperbolic arc tangent of the argument
 */
template <typename T>
inline stan::math::complex<fvar<T>> atanh(
    const stan::math::complex<fvar<T>>& z) {
  return internal::complex_atanh(z);
}

}  // namespace math
}  // namespace stan
#endif
