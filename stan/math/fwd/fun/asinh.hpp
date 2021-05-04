#ifndef STAN_MATH_FWD_FUN_ASINH_HPP
#define STAN_MATH_FWD_FUN_ASINH_HPP

#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/meta.hpp>
#include <stan/math/prim/fun/asinh.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {

template <typename T>
inline fvar<T> asinh(const fvar<T>& x) {
  using std::sqrt;
  return fvar<T>(asinh(x.val_), x.d_ / sqrt(square(x.val_) + 1));
}

/**
 * Return the hyperbolic arcsine of the complex argument.
 *
 * @tparam T autodiff value type
 * @param[in] z argument
 * @return hyperbolic arcsine of the argument
 */
template <typename T>
inline std::complex<fvar<T>> asinh(const std::complex<fvar<T>>& z) {
  return internal::complex_asinh(z);
}

}  // namespace math
}  // namespace stan
#endif
