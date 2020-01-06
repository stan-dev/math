#ifndef STAN_MATH_PRIM_MAT_FUN_COS_HPP
#define STAN_MATH_PRIM_MAT_FUN_COS_HPP

#include <stan/math/prim/vectorize/apply_scalar_unary.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Structure to wrap cos() so it can be vectorized.
 *
 * @tparam T type of variable
 * @param x angle in radians
 * @return Cosine of x.
 */
struct cos_fun {
  template <typename T>
  static inline T fun(const T& x) {
    using std::cos;
    return cos(x);
  }
};

/**
 * Vectorized version of cos().
 *
 * @tparam T type of container
 * @param x angles in radians
 * @return Cosine of each value in x.
 */
template <typename T>
inline typename apply_scalar_unary<cos_fun, T>::return_t cos(const T& x) {
  return apply_scalar_unary<cos_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
