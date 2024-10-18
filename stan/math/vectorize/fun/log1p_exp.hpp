
#ifndef STAN_MATH_VECTORIZED_FUN_LOG1P_EXP_HPP
#define STAN_MATH_VECTORIZED_FUN_LOG1P_EXP_HPP
#include <stan/math/prim/fun/log1p_exp.hpp>
namespace stan {
namespace math {

/**
 * Vectorized version of log1p_exp().
 *
 * @tparam T type of container
 * @param x container
 * @return Natural log of (1 + exp()) applied to each value in x.
 */
template <typename T,
          require_not_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr,
          require_not_var_matrix_t<T>* = nullptr>
inline auto log1p_exp(const T& x) {
  return apply_scalar_unary<log1p_exp_fun, T>::apply(x);
}


} // namespace math
} // namespace stan
#endif 

