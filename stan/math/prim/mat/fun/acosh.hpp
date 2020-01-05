#ifndef STAN_MATH_PRIM_MAT_FUN_ACOSH_HPP
#define STAN_MATH_PRIM_MAT_FUN_ACOSH_HPP

#include <stan/math/prim/vectorize/apply_scalar_unary.hpp>
#include <stan/math/prim/scal/fun/acosh.hpp>

namespace stan {
namespace math {

/**
 * Structure to wrap acosh() so it can be vectorized.
 */
struct acosh_fun {
  /**
   * Return the inverse hypberbolic cosine of the specified argument.
   *
   * @tparam T type of argument
   * @param x argument
   * @return Inverse hyperbolic cosine of the argument.
   */
  template <typename T>
  static inline T fun(const T& x) {
    return acosh(x);
  }
};

/**
 * Return the elementwise application of <code>acosh()</code> to
 * specified argument container.  The return type promotes the
 * underlying scalar argument type to double if it is an integer,
 * and otherwise is the argument type.
 *
 * @tparam T type of container
 * @param x container
 * @return Elementwise acosh of members of container.
 */
template <typename T>
inline typename apply_scalar_unary<acosh_fun, T>::return_t acosh(const T& x) {
  return apply_scalar_unary<acosh_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
