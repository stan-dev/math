#ifndef STAN_MATH_PRIM_FUN_ASINH_HPP
#define STAN_MATH_PRIM_FUN_ASINH_HPP

#include <stan/math/prim/meta.hpp>
#include <cmath>

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

}  // namespace math
}  // namespace stan

#endif
