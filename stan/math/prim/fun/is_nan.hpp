#ifndef STAN_MATH_PRIM_FUN_IS_NAN_HPP
#define STAN_MATH_PRIM_FUN_IS_NAN_HPP

#include <stan/math/prim/meta.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Returns true if the input is NaN and false otherwise.
 *
 * Delegates to <code>std::isnan</code>.
 *
 * @param x Value to test.
 * @return <code>true</code> if the value is NaN.
 */
template <typename T, typename = require_arithmetic_t<T>>
inline bool is_nan(T x) {
  return std::isnan(x);
}

template <typename T, typename = require_eigen_t<T>>
inline bool is_nan(const T& x) {
  return x.hasNan();
}

}  // namespace math
}  // namespace stan

#endif
