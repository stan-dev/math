#ifndef STAN_MATH_PRIM_MAT_FUN_LOG1M_EXP_HPP
#define STAN_MATH_PRIM_MAT_FUN_LOG1M_EXP_HPP

#include <stan/math/prim/vectorize/apply_scalar_unary.hpp>
#include <stan/math/prim/scal/fun/log1m_exp.hpp>

namespace stan {
namespace math {

/**
 * Structure to wrap log1m_exp() so it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return Natural log of (1 - exp(x)).
 */
struct log1m_exp_fun {
  template <typename T>
  static inline T fun(const T& x) {
    return log1m_exp(x);
  }
};

/**
 * Vectorized version of log1m_exp().
 *
 * @tparam T type of container
 * @param x container
 * @return Natural log of (1 - exp()) applied to each value in x.
 */
template <typename T>
inline typename apply_scalar_unary<log1m_exp_fun, T>::return_t log1m_exp(
    const T& x) {
  return apply_scalar_unary<log1m_exp_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
