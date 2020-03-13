#ifndef STAN_MATH_PRIM_FUN_CBRT_HPP
#define STAN_MATH_PRIM_FUN_CBRT_HPP

#include <stan/math/prim/meta.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the cube root of the specified value
 *
 * @param[in] x Argument.
 * @return Cube root of the argument.
 * @throw std::domain_error If argument is negative.
 */
inline double cbrt(double x) { return std::cbrt(x); }

/**
 * Integer version of cbrt.
 *
 * @param[in] x Argument.
 * @return Cube root of the argument.
 * @throw std::domain_error If argument is less than 1.
 */
inline double cbrt(int x) { return std::cbrt(x); }

/**
 * Structure to wrap cbrt() so it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return Cube root of x.
 */
struct cbrt_fun {
  template <typename T>
  static inline T fun(const T& x) {
    return cbrt(x);
  }
};

/**
 * Vectorized version of cbrt().
 *
 * @tparam T type of container
 * @param x container
 * @return Cube root of each value in x.
 */
template <typename T>
inline auto cbrt(const T& x) {
  return apply_scalar_unary<cbrt_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
