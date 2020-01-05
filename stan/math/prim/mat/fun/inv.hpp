#ifndef STAN_MATH_PRIM_MAT_FUN_INV_HPP
#define STAN_MATH_PRIM_MAT_FUN_INV_HPP

#include <stan/math/prim/vectorize/apply_scalar_unary.hpp>
#include <stan/math/prim/scal/fun/inv.hpp>

namespace stan {
namespace math {

/**
 * Structure to wrap inv() so that it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return 1 / x.
 */
struct inv_fun {
  template <typename T>
  static inline T fun(const T& x) {
    return inv(x);
  }
};

/**
 * Vectorized version of inv().
 *
 * @tparam T type of container
 * @param x container
 * @return 1 divided by each value in x.
 */
template <typename T>
inline typename apply_scalar_unary<inv_fun, T>::return_t inv(const T& x) {
  return apply_scalar_unary<inv_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
