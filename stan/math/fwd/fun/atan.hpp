#ifndef STAN_MATH_FWD_FUN_ATAN_HPP
#define STAN_MATH_FWD_FUN_ATAN_HPP

#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/meta.hpp>
#include <stan/math/prim/fun/atan.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <cmath>
#include <stan/math/prim/core/complex_base.hpp>

namespace stan {
namespace math {

template <typename T>
inline fvar<T> atan(const fvar<T>& x) {
  using std::atan;
  return fvar<T>(atan(x.val_), x.d_ / (1 + square(x.val_)));
}

/**
 * Return the arc tangent of the complex argument.
 *
 * @tparam T autodiff value type
 * @param[in] z argument
 * @return arc tanget of the argument
 */
template <typename T>
inline stan::math::complex<fvar<T>> atan(
    const stan::math::complex<fvar<T>>& z) {
  return internal::complex_atan(z);
}

}  // namespace math
}  // namespace stan
#endif
