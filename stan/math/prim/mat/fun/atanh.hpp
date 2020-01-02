#ifndef STAN_MATH_PRIM_MAT_FUN_ATANH_HPP
#define STAN_MATH_PRIM_MAT_FUN_ATANH_HPP

#include <stan/math/prim/vectorize/apply_scalar_unary.hpp>
#include <stan/math/prim/scal/fun/atanh.hpp>

namespace stan {
namespace math {

/**
 * Structure to wrap atanh() so it can be vectorized.
 */
struct atanh_fun {
  /**
   * Return the inverse hypberbolic tangent of the specified argument.
   *
   * @tparam T type of argument
   * @param x argument
   * @return Inverse hyperbolic tangent of the argument.
   */
  template <typename T>
  static inline T fun(const T& x) {
    return atanh(x);
  }
};

/**
 * Return the elementwise application of <code>atanh()</code> to
 * specified argument container.  The return type promotes the
 * underlying scalar argument type to double if it is an integer,
 * and otherwise is the argument type.
 *
 * @tparam T type of container
 * @param x container
 * @return Elementwise atanh of members of container.
 */
template <typename T>
inline typename apply_scalar_unary<atanh_fun, T>::return_t atanh(const T& x) {
  return apply_scalar_unary<atanh_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
