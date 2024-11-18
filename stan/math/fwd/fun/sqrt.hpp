#ifndef STAN_MATH_FWD_FUN_SQRT_HPP
#define STAN_MATH_FWD_FUN_SQRT_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/cos.hpp>
#include <stan/math/fwd/fun/inv_sqrt.hpp>
#include <stan/math/fwd/fun/sin.hpp>
#include <stan/math/fwd/fun/hypot.hpp>
#include <stan/math/prim/fun/sqrt.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {

template <typename T>
inline fvar<T> sqrt(const fvar<T>& x) {
  if (value_of_rec(x.val_) == 0.0) {
    return fvar<T>(sqrt(x.val_), 0.0 * x.d_);
  }
  return fvar<T>(sqrt(x.val_), 0.5 * x.d_ * inv_sqrt(x.val_));
}

/**
 * Return the square root of the complex argument.
 *
 * @tparam T autodiff value type
 * @param[in] z argument
 * @return square root of the argument
 */
template <typename T>
inline std::complex<fvar<T>> sqrt(const std::complex<fvar<T>>& z) {
  return internal::complex_sqrt(z);
}

}  // namespace math
}  // namespace stan
#endif
