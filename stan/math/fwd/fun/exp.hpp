#ifndef STAN_MATH_FWD_FUN_EXP_HPP
#define STAN_MATH_FWD_FUN_EXP_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {
template <typename T>
inline fvar<T> exp(const fvar<T>& x) {
  using std::exp;
  return fvar<T>(exp(x.val_), x.d_ * exp(x.val_));
}

// template <typename T>
// inline std::complex<fvar<T>> exp(const std::complex<fvar<T>>& z) {
//  return complex_exp(z);
// }

}  // namespace math
}  // namespace stan
#endif
