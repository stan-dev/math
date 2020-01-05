#ifndef STAN_MATH_PRIM_MAT_FUN_INV_LOGIT_HPP
#define STAN_MATH_PRIM_MAT_FUN_INV_LOGIT_HPP

#include <stan/math/prim/vectorize/apply_scalar_unary.hpp>
#include <stan/math/prim/scal/fun/inv_logit.hpp>

namespace stan {
namespace math {

/**
 * Structure to wrap inv_logit() so that it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return Inverse logit of x.
 */
struct inv_logit_fun {
  template <typename T>
  static inline T fun(const T& x) {
    return inv_logit(x);
  }
};

/**
 * Vectorized version of inv_logit().
 *
 * @tparam T type of container
 * @param x container
 * @return Inverse logit applied to each value in x.
 */
template <typename T>
inline typename apply_scalar_unary<inv_logit_fun, T>::return_t inv_logit(
    const T& x) {
  return apply_scalar_unary<inv_logit_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
