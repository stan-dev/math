#ifndef STAN_MATH_FWD_FUN_POLAR_HPP
#define STAN_MATH_FWD_FUN_POLAR_HPP

#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/meta.hpp>
#include <stan/math/prim/fun/polar.hpp>
#include <complex>
#include <type_traits>

namespace stan {
namespace math {

/**
 * Returns complex number with specified magnitude and phase angle.
 *
 * @tparam T autodiff value type
 * @param[in] r magnitude
 * @param[in] theta phase angle
 * @return complex number with magnitude and phase angle
 */
template <typename T>
inline std::complex<fvar<T>> polar(const fvar<T>& r, const fvar<T>& theta) {
  return internal::complex_polar(r, theta);
}

/**
 * Returns complex number with specified magnitude and phase angle.
 *
 * @tparam T autodiff value type for magnitude
 * @tparam U arithmetic type for phase angle
 * @param[in] r magnitude
 * @param[in] theta phase angle
 * @return complex number with magnitude and phase angle
 */
template <typename T, typename U>
inline std::complex<fvar<T>> polar(const fvar<T>& r, U theta) {
  return internal::complex_polar(r, theta);
}

/**
 * Returns complex number with specified magnitude and phase angle.
 *
 * @tparam T autodiff value type for phase angle
+* * @tparam U arithmetic type for magnitude
 * @param[in] r magnitude
 * @param[in] theta phase angle
 * @return complex number with magnitude and phase angle
 */
template <typename T, typename U>
inline std::complex<fvar<T>> polar(U r, const fvar<T>& theta) {
  return internal::complex_polar(r, theta);
}
}  // namespace math
}  // namespace stan

#endif
