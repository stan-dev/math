#ifndef STAN_MATH_PRIM_FUN_IS_NAN_HPP
#define STAN_MATH_PRIM_FUN_IS_NAN_HPP

#include <stan/math/prim/meta.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Tests if the input is Not a Number (NaN)
 *
 * Integer specialization, always return false.
 *
 * @tparam T Integer type.
 * @param x Value to test.
 * @return <code>false</code>.
 */
template <typename T, require_integral_t<T>* = nullptr>
inline bool is_nan(T x) {
  return false;
}

/**
 * Tests if the input is Not a Number (NaN)
 *
 * Delegates to <code>std::isnan</code>.
 *
 * @tparam T Floating point type.
 * @param x Value to test.
 * @return <code>true</code> if the value is NaN.
 */
template <typename T, require_floating_point_t<T>* = nullptr>
inline bool is_nan(T x) {
  return std::isnan(x);
}

template <typename T, require_eigen_t<T>* = nullptr>
inline bool is_nan(const T& x) {
  return x.hasNan();
}

}  // namespace math
}  // namespace stan

#endif
