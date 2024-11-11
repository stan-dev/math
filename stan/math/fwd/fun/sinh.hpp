#ifndef STAN_MATH_FWD_FUN_SINH_HPP
#define STAN_MATH_FWD_FUN_SINH_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/cosh.hpp>
#include <stan/math/fwd/fun/exp.hpp>
#include <stan/math/prim/fun/sinh.hpp>
#include <complex>

namespace stan {
namespace math {

template <typename T>
inline fvar<T> sinh(const fvar<T>& x) {
  return fvar<T>(sinh(x.val_), x.d_ * cosh(x.val_));
}

/**
 * Return the hyperbolic sine of the complex argument.
 *
 * @tparam T autodiff value type
 * @param[in] z argument
 * @return hyperbolic sine of the argument
 */
template <typename T>
inline std::complex<fvar<T>> sinh(const std::complex<fvar<T>>& z) {
  return internal::complex_sinh(z);
}

}  // namespace math
}  // namespace stan
#endif
