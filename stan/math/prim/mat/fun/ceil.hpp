#ifndef STAN_MATH_PRIM_MAT_FUN_CEIL_HPP
#define STAN_MATH_PRIM_MAT_FUN_CEIL_HPP

#include <stan/math/prim/vectorize/apply_scalar_unary.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Structure to wrap ceil() so it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return Least integer >= x.
 */
struct ceil_fun {
  template <typename T>
  static inline T fun(const T& x) {
    using std::ceil;
    return ceil(x);
  }
};

/**
 * Vectorized version of ceil().
 *
 * @tparam T type of container
 * @param x container
 * @return Least integer >= each value in x.
 */
template <typename T>
inline typename apply_scalar_unary<ceil_fun, T>::return_t ceil(const T& x) {
  return apply_scalar_unary<ceil_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
