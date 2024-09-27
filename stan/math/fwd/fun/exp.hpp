#ifndef STAN_MATH_FWD_FUN_EXP_HPP
#define STAN_MATH_FWD_FUN_EXP_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <cmath>
#include <stan/math/prim/core/complex_base.hpp>

namespace stan {
namespace math {
template <typename T>
inline fvar<T> exp(const fvar<T>& x) {
  using std::exp;
  return fvar<T>(exp(x.val_), x.d_ * exp(x.val_));
}

/**
 * Return the natural exponentiation (base e) of the specified complex number.
 *
 * @tparam T value type of autodiff variable
 * @param z complex argument
 * @return exponentiation of argument
 */
template <typename T>
inline stan::math::complex<fvar<T>> exp(const stan::math::complex<fvar<T>>& z) {
  return internal::complex_exp(z);
}

}  // namespace math
}  // namespace stan
#endif
