#ifndef STAN_MATH_PRIM_MAT_FUN_SIN_HPP
#define STAN_MATH_PRIM_MAT_FUN_SIN_HPP

#include <stan/math/prim/vectorize/apply_scalar_unary.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Structure to wrap sin() so it can be vectorized.
 *
 * @tparam T type of argument
 * @param x angle in radians
 * @return Sine of x.
 */
struct sin_fun {
  template <typename T>
  static inline T fun(const T& x) {
    using std::sin;
    return sin(x);
  }
};

/**
 * Vectorized version of sin().
 *
 * @tparam T type of container
 * @param x angles in radians
 * @return Sine of each value in x.
 */
template <typename T>
inline typename apply_scalar_unary<sin_fun, T>::return_t sin(const T& x) {
  return apply_scalar_unary<sin_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
