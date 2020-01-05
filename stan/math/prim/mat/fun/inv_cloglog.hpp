#ifndef STAN_MATH_PRIM_MAT_FUN_INV_CLOGLOG_HPP
#define STAN_MATH_PRIM_MAT_FUN_INV_CLOGLOG_HPP

#include <stan/math/prim/vectorize/apply_scalar_unary.hpp>
#include <stan/math/prim/scal/fun/inv_cloglog.hpp>

namespace stan {
namespace math {

/**
 * Structure to wrap inv_cloglog() so that it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return 1 - exp(-exp(x)).
 */
struct inv_cloglog_fun {
  template <typename T>
  static inline T fun(const T& x) {
    return inv_cloglog(x);
  }
};

/**
 * Vectorized version of inv_cloglog().
 *
 * @tparam T type of container
 * @param x container
 * @return 1 - exp(-exp()) applied to each value in x.
 */
template <typename T>
inline typename apply_scalar_unary<inv_cloglog_fun, T>::return_t inv_cloglog(
    const T& x) {
  return apply_scalar_unary<inv_cloglog_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
