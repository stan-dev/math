#ifndef STAN_MATH_FWD_FUN_COS_HPP
#define STAN_MATH_FWD_FUN_COS_HPP

#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/meta.hpp>
#include <stan/math/prim/fun/cos.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {

template <typename T>
inline fvar<T> cos(const fvar<T>& x) {
  using std::cos;
  using std::sin;
  return fvar<T>(cos(x.val_), x.d_ * -sin(x.val_));
}

/**
 * Return the cosine of the complex argument.
 *
 * @tparam T autodiff value type
 * @param[in] z argument
 * @return cosine of the argument
 */
template <typename T>
inline std::complex<fvar<T>> cos(const std::complex<fvar<T>>& z) {
  return internal::complex_cos(z);
}

}  // namespace math
}  // namespace stan
#endif
