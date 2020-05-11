#ifndef STAN_MATH_PRIM_FUN_FMA_HPP
#define STAN_MATH_PRIM_FUN_FMA_HPP

#include <stan/math/prim/meta.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the product of the first two arguments plus the third
 * argument.
 *
 * <p><i>Warning:</i> This does not delegate to the high-precision
 * platform-specific <code>fma()</code> implementation.
 *
 * @param x First argument.
 * @param y Second argument.
 * @param z Third argument.
 * @return The product of the first two arguments plus the third
 * argument.
 */
template <typename T1, typename T2, typename T3,
          typename = require_all_arithmetic_t<T1, T2, T3>>
inline double fma(T1 x, T2 y, T3 z) {
  using std::fma;
  return fma(x, y, z);
}

template <typename T1, typename T2, typename T3,
 require_any_eigen_t<T1, T2, T3>* = nullptr,
 require_all_vt_arithmetic<T1, T2, T3>* = nullptr>
inline auto fma(const T1& x, const T2& y, const T3& z) {
  using std::fma;
  return (as_array_or_scalar(x) * as_array_or_scalar(y) + as_array_or_scalar(z)).matrix().eval();
}


}  // namespace math
}  // namespace stan
#endif
