#ifndef STAN_MATH_PRIM_FUN_POLAR_HPP
#define STAN_MATH_PRIM_FUN_POLAR_HPP

#include <stan/math/prim/fun/cos.hpp>
#include <stan/math/prim/fun/is_inf.hpp>
#include <stan/math/prim/meta.hpp>
#include <cmath>
#include <complex>
#include <limits>

namespace stan {
namespace math {
namespace internal {
/**
 * Returns complex number with specified magnitude and phase angle.
 *
 * @param[in] r magnitude
 * @param[in] theta phase angle
 * @return complex number with magnitude and phase angle
 */
template <typename U, typename V>
inline complex_return_t<U, V> complex_polar(const U& r, const V& theta) {
  using std::cos;
  using std::sin;
  if (!(r >= 0) || is_inf(theta)) {
    return {std::numeric_limits<double>::quiet_NaN()};
  }
  return {r * cos(theta), r * sin(theta)};
}
}  // namespace internal

/**
 * Returns the complex number with specified magnitude and phase angle.
 *
 * @param[in] r magnitude
 * @param[in] theta phase angle
 * @return complex number with magnitude and phase angle
 */
inline std::complex<double> polar(double r, double theta) {
  return internal::complex_polar(r, theta);
}
inline std::complex<double> polar(double r, int theta) {
  return internal::complex_polar(r, static_cast<double>(theta));
}
inline std::complex<double> polar(int r, double theta) {
  return internal::complex_polar(static_cast<double>(r), theta);
}
inline std::complex<double> polar(int r, int theta) {
  return internal::complex_polar(static_cast<double>(r),
                                 static_cast<double>(theta));
}

}  // namespace math
}  // namespace stan

#endif
