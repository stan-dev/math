#ifndef STAN_MATH_PRIM_MAT_FUN_INV_SQRT_HPP
#define STAN_MATH_PRIM_MAT_FUN_INV_SQRT_HPP

#include <stan/math/prim/vectorize/apply_scalar_unary.hpp>
#include <stan/math/prim/scal/fun/inv_sqrt.hpp>

namespace stan {
namespace math {

/**
 * Structure to wrap inv_sqrt() so that it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return 1 / sqrt of x.
 */
struct inv_sqrt_fun {
  template <typename T>
  static inline T fun(const T& x) {
    return inv_sqrt(x);
  }
};

/**
 * Vectorized version of inv_sqrt().
 *
 * @tparam T type of container
 * @param x container
 * @return 1 / sqrt of each value in x.
 */
template <typename T>
inline typename apply_scalar_unary<inv_sqrt_fun, T>::return_t inv_sqrt(
    const T& x) {
  return apply_scalar_unary<inv_sqrt_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
