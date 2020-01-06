#ifndef STAN_MATH_PRIM_MAT_FUN_ACOS_HPP
#define STAN_MATH_PRIM_MAT_FUN_ACOS_HPP

#include <stan/math/prim/vectorize/apply_scalar_unary.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Structure to wrap acos() so it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return Arc cosine of variable in radians.
 */
struct acos_fun {
  template <typename T>
  static inline T fun(const T& x) {
    using std::acos;
    return acos(x);
  }
};

/**
 * Vectorized version of acos().
 *
 * @tparam T type of container
 * @param x container
 * @return Arc cosine of each variable in the container, in radians.
 */
template <typename T>
inline typename apply_scalar_unary<acos_fun, T>::return_t acos(const T& x) {
  return apply_scalar_unary<acos_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
