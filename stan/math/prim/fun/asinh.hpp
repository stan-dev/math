#ifndef STAN_MATH_PRIM_FUN_ASINH_HPP
#define STAN_MATH_PRIM_FUN_ASINH_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/copysign.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {

/**
 * Return the inverse hyperbolic sine of the specified value.
 * Returns infinity for infinity argument and -infinity for
 * -infinity argument.
 * Returns nan for nan argument.
 *
 * @param[in] x Argument.
 * @return Inverse hyperbolic sine of the argument.
 */
inline double asinh(double x) { return std::asinh(x); }

/**
 * Integer version of asinh.
 *
 * @param[in] x Argument.
 * @return Inverse hyperbolic sine of the argument.
 */
inline double asinh(int x) { return std::asinh(x); }

/**
 * Structure to wrap asinh() so it can be vectorized.
 *
 * @tparam T argument scalar type
 * @param x argument
 * @return inverse hyperbolic sine of argument in radians.
 */
struct asinh_fun {
  template <typename T>
  static inline T fun(const T& x) {
    return asinh(x);
  }
};

/**
 * Vectorized version of asinh().
 *
 * @tparam T type of container
 * @param x container
 * @return Inverse hyperbolic sine of each value in the container.
 */
template <typename T>
inline auto asinh(const T& x) {
  return apply_scalar_unary<asinh_fun, T>::apply(x);
}

// namespace internal {
// /**
//  * Return the hyperbolic arc sine of the complex argument.
//  *
//  * @tparam V value type of argument
//  * @param[in] z argument
//  * @return hyperbolic arc sine of the argument
//  */
// template <typename V>
// inline std::complex<V> complex_asinh(const std::complex<V>& z) {
//   std::complex<double> y_d = asinh(value_of_rec(z));
//   auto y = log(z + sqrt(1 + z * z));
//   return copysign(y, y_d);
// }
// }  // namespace internal

}  // namespace math
}  // namespace stan

#endif
