
#ifndef STAN_MATH_VECTORIZED_FUN_LOG1M_INV_LOGIT_HPP
#define STAN_MATH_VECTORIZED_FUN_LOG1M_INV_LOGIT_HPP
#include <stan/math/prim/fun/log1m_inv_logit.hpp>
namespace stan {
namespace math {

/**
 * Return the elementwise application of
 * <code>log1m_inv_logit()</code> to specified argument container.
 * The return type promotes the underlying scalar argument type to
 * double if it is an integer, and otherwise is the argument type.
 *
 * @tparam T type of container
 * @param x container
 * @return Elementwise log1m_inv_logit of members of container.
 */
template <typename T, require_not_var_matrix_t<T>* = nullptr,
          require_not_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr>
inline typename apply_scalar_unary<log1m_inv_logit_fun, T>::return_t
log1m_inv_logit(const T& x) {
  return apply_scalar_unary<log1m_inv_logit_fun, T>::apply(x);
}


} // namespace math
} // namespace stan
#endif 

