#ifndef STAN_MATH_REV_FUN_POLAR_HPP
#define STAN_MATH_REV_FUN_POLAR_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/is_inf.hpp>
#include <stan/math/rev/fun/cos.hpp>
#include <stan/math/rev/fun/sin.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/isinf.hpp>
#include <stan/math/prim/fun/isfinite.hpp>
#include <stan/math/prim/fun/is_inf.hpp>
#include <stan/math/prim/fun/polar.hpp>
#include <complex>
#include <type_traits>

namespace stan {
namespace math {

/**
 * Returns complex number with specified magnitude and phase angle.
 *
 * @param[in] r magnitude
 * @param[in] theta phase angle
 * @return complex number with magnitude and phase angle
 */
inline std::complex<var> polar(const var& r, const var& theta) {
  return internal::complex_polar(r, theta);
}

/**
 * Returns complex number with specified magnitude and phase angle.
 *
 * @tparam T arithmetic type of magnitude
 * @param[in] r magnitude
 * @param[in] theta phase angle
 * @return complex number with magnitude and phase angle
 */
template <typename T>
inline std::complex<var> polar(T r, const var& theta) {
  return internal::complex_polar(r, theta);
}

/**
 * Returns complex number with specified magnitude and phase angle.
 *
 * @tparam T arithmetic type of phase angle
 * @param[in] r magnitude
 * @param[in] theta phase angle
 * @return complex number with magnitude and phase angle
 */
template <typename T>
inline std::complex<var> polar(const var& r, T theta) {
  return internal::complex_polar(r, theta);
}

}  // namespace math
}  // namespace stan

#endif
