#ifndef STAN_MATH_PRIM_FUN_POW_HPP
#define STAN_MATH_PRIM_FUN_POW_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {

// This overload is required becuase MINGW32 does not include the
// specialization std::pow(double, int) as required by the C++11
// standard library interface specification.  Use with other compilers
// than MINGW32 would introduce an ambiguity with the standard library.
#if __MINGW32__
/**
 * Return the first argument raised to the power of the second
 * argument.
 *
 * @param x base
 * @param y exponent
 * @return base raised to the power of the exponent
 */
double pow(double x, int y) { return std::pow(x, static_cast<double>(y)); }
#endif

namespace internal {

/**
 * Return the first argument raised to the power of the second
 * argument.  At least one of the arguments must be a complex number.
 *
 * @tparam U type of base
 * @tparam V type of exponent
 * @param[in] x base
 * @param[in] y exponent
 * @return base raised to the power of the exponent
 */
template <typename U, typename V>
inline complex_return_t<U, V> complex_pow(const U& x, const V& y) {
  return exp(y * log(x));
}
}  // namespace internal
}  // namespace math
}  // namespace stan
#endif
