#ifndef STAN_MATH_FWD_FUN_COSH_HPP
#define STAN_MATH_FWD_FUN_COSH_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/exp.hpp>
#include <stan/math/prim/fun/cosh.hpp>
#include <cmath>
#include <stan/math/prim/core/complex_base.hpp>

namespace stan {
namespace math {

template <typename T>
inline fvar<T> cosh(const fvar<T>& x) {
  using std::cosh;
  using std::sinh;
  return fvar<T>(cosh(x.val_), x.d_ * sinh(x.val_));
}

/**
 * Return the hyperbolic cosine of the complex argument.
 *
 * @tparam T autodiff value type
 * @param[in] z argument
 * @return hyperbolic cosine of the argument
 */
template <typename T>
inline stan::math::complex<fvar<T>> cosh(
    const stan::math::complex<fvar<T>>& z) {
  return stan::math::internal::complex_cosh(z);
}

}  // namespace math
}  // namespace stan
#endif
