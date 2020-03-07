#ifndef STAN_MATH_PRIM_FUN_EXP2_HPP
#define STAN_MATH_PRIM_FUN_EXP2_HPP

#include <stan/math/prim/meta.hpp>
#include <cmath>

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
template <typename T, typename = require_arithmetic_t<T>>
inline double exp2(T y) {
  using std::exp2;
  return exp2(y);
}

/**
 * Structure to wrap exp2() so it can be vectorized.
 */
struct exp2_fun {
  /**
   * Return the base two exponent of the specified argument.
   *
   * @tparam T type of argument
   * @param x argument
   * @return Base two exponent of the argument.
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
 * @tparam T type of container
 * @param x container
 * @return Elementwise exp2 of members of container.
 */
template <typename T, typename = require_vector_like_t<T>>
inline auto exp2(const T& x) {
  return apply_scalar_unary<exp2_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
