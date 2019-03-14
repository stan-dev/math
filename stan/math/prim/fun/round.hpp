#ifndef STAN_MATH_PRIM_FUN_ROUND_HPP
#define STAN_MATH_PRIM_FUN_ROUND_HPP

#include <cmath>
#include <stan/math/prim/vectorize/apply_scalar_unary.hpp>




namespace stan {
namespace math {

/**
 * Return the closest integer to the specified argument, with
 * halfway cases rounded away from zero.
 *
 * @param x Argument.
 * @return The rounded value of the argument.
 */
inline double round(double x) { return std::round(x); }

/**
 * Return the closest integer to the specified argument, with
 * halfway cases rounded away from zero.
 *
 * @param x Argument.
 * @return The rounded value of the argument.
 */
inline double round(int x) { return std::round(x); }













/**
 * Structure to wrap round() so it can be vectorized.
 * @param x Argument variable.
 * @tparam T Argument type.
 * @return Rounded value of x.
 */
struct round_fun {
  template <typename T>
  static inline T fun(const T& x) {
    using stan::math::round;
    return round(x);
  }
};

/**
 * Vectorized version of round.
 * @param x Container.
 * @tparam T Container type.
 * @return Rounded value of each value in x.
 */
template <typename T>
inline typename apply_scalar_unary<round_fun, T>::return_t round(const T& x) {
  return apply_scalar_unary<round_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
