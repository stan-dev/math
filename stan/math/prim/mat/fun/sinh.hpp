#ifndef STAN_MATH_PRIM_MAT_FUN_SINH_HPP
#define STAN_MATH_PRIM_MAT_FUN_SINH_HPP

#include <stan/math/prim/vectorize/apply_scalar_unary.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Structure to wrap sinh() so that it can be vectorized.
 *
 * @tparam T type of argument
 * @param x angle in radians
 * @return Hyperbolic sine of x.
 */
struct sinh_fun {
  template <typename T>
  static inline T fun(const T& x) {
    using std::sinh;
    return sinh(x);
  }
};

/**
 * Vectorized version of sinh().
 *
 * @tparam T type of container
 * @param x container
 * @return Hyperbolic sine of each variable in x.
 */
template <typename T>
inline typename apply_scalar_unary<sinh_fun, T>::return_t sinh(const T& x) {
  return apply_scalar_unary<sinh_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
