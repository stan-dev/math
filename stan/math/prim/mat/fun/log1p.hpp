#ifndef STAN_MATH_PRIM_MAT_FUN_LOG1P_HPP
#define STAN_MATH_PRIM_MAT_FUN_LOG1P_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/fun/log1p.hpp>
#include <stan/math/prim/vectorize/apply_scalar_unary.hpp>

namespace stan {
namespace math {

/**
 * Structure to wrap log1p() so it can be vectorized.
 */
struct log1p_fun {
  /**
   * Return the natural logarithm of one plus the specified value.
   *
   * @tparam T type of argument
   * @param x argument
   * @return natural log of one plus the argument
   */
  template <typename T>
  static inline T fun(const T& x) {
    return log1p(x);
  }
};

/**
 * Return the elementwise application of <code>log1p()</code> to
 * specified argument container.  The return type promotes the
 * underlying scalar argument type to double if it is an integer,
 * and otherwise is the argument type.
 *
 * @tparam T type of container
 * @param x container
 * @return Elementwise log1p of members of container.
 */
template <typename T>
inline auto log1p(const T& x) {
  return apply_scalar_unary<log1p_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
