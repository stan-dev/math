#ifndef STAN_MATH_FWD_FUN_ACOSH_HPP
#define STAN_MATH_FWD_FUN_ACOSH_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/fun/acosh.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/fun/sqrt.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {

template <typename T>
inline fvar<T> acosh(const fvar<T>& x) {
  using std::sqrt;
  return fvar<T>(acosh(x.val_), x.d_ / sqrt(square(x.val_) - 1));
}

/**
 * Return the hyperbolic arc cosine of the complex argument.
 *
 * @tparam T autodiff value type
 * @param[in] z argument
 * @return hyperbolic arc cosine of the argument
 */
template <typename T>
inline std::complex<fvar<T>> acosh(const std::complex<fvar<T>>& z) {
  return internal::complex_acosh(z);
}

}  // namespace math
}  // namespace stan
#endif
