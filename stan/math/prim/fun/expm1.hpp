#ifndef STAN_MATH_PRIM_FUN_EXPM1_HPP
#define STAN_MATH_PRIM_FUN_EXPM1_HPP

#include <stan/math/prim/meta.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the natural exponentiation of x minus one.
 * Returns infinity for infinity argument and -infinity for
 * -infinity argument.
 *
 * @param[in] x Argument.
 * @return Natural exponentiation of argument minus one.
 */
inline double expm1(double x) { return std::expm1(x); }

/**
 * Integer version of expm1.
 *
 * @param[in] x Argument.
 * @return Natural exponentiation of argument minus one.
 */
inline double expm1(int x) { return std::expm1(x); }

/**
 * Structure to wrap expm1() so that it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return Natural exponential of x minus one.
 */
struct expm1_fun {
  template <typename T>
  static inline T fun(const T& x) {
    return expm1(x);
  }
};

/**
 * Vectorized version of expm1().
 *
 * @tparam T type of container
 * @param x container
 * @return Natural exponential of each value in x minus one.
 */
template <typename T>
inline auto expm1(const T& x) {
  return apply_scalar_unary<expm1_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
