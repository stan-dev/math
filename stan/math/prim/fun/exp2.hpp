#ifndef STAN_MATH_PRIM_FUN_EXP2_HPP
#define STAN_MATH_PRIM_FUN_EXP2_HPP

#include <boost/math/tools/promotion.hpp>
#include <cmath>
#include <stan/math/prim/vectorize/apply_scalar_unary.hpp>






namespace stan {
namespace math {

/**
 * Return the exponent base 2 of the specified argument (C99,
 * C++11).
 *
 * The exponent base 2 function is defined by
 *
 * <code>exp2(y) = pow(2.0, y)</code>.
 *
 * @param y argument.
 * @return exponent base 2 of argument.
 */
inline double exp2(double y) {
  using std::pow;
  return pow(2.0, y);
}

/**
 * Return the exponent base 2 of the specified argument (C99,
 * C++11).
 *
 * @param y argument
 * @return exponent base 2 of argument
 */
inline double exp2(int y) { return exp2(static_cast<double>(y)); }













/**
 * Structure to wrap exp2() so it can be vectorized.
 */
struct exp2_fun {
  /**
   * Return the base two exponent of the specified argument.
   *
   * @param x Argument.
   * @return Base two exponent of the argument.
   * @tparam T Argument type.
   */
  template <typename T>
  static inline T fun(const T& x) {
    return exp2(x);
  }
};

/**
 * Return the elementwise application of <code>exp2()</code> to
 * specified argument container.  The return type promotes the
 * underlying scalar argument type to double if it is an integer,
 * and otherwise is the argument type.
 *
 * @tparam T Container type.
 * @param x Container.
 * @return Elementwise exp2 of members of container.
 */
template <typename T>
inline typename apply_scalar_unary<exp2_fun, T>::return_t exp2(const T& x) {
  return apply_scalar_unary<exp2_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
