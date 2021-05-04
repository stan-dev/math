#ifndef STAN_MATH_FWD_FUN_ACOS_HPP
#define STAN_MATH_FWD_FUN_ACOS_HPP

#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/meta.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/fun/acos.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {

template <typename T>
inline fvar<T> acos(const fvar<T>& x) {
  using std::acos;
  using std::sqrt;
  return fvar<T>(acos(x.val_), x.d_ / -sqrt(1 - square(x.val_)));
}

/**
 * Return the arc cosine of the complex argument.
 *
 * @tparam T autodiff value type
 * @param x argument
 * @return arc cosine of the argument
 */
template <typename T>
inline std::complex<fvar<T>> acos(const std::complex<fvar<T>>& x) {
  return internal::complex_acos(x);
}

}  // namespace math
}  // namespace stan
#endif
