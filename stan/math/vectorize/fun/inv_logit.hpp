
#ifndef STAN_MATH_VECTORIZED_FUN_INV_LOGIT_HPP
#define STAN_MATH_VECTORIZED_FUN_INV_LOGIT_HPP
#include <stan/math/prim/fun/inv_logit.hpp>
namespace stan {
namespace math {

/**
 * Vectorized version of inv_logit().
 *
 * @tparam T type of container
 * @param x container
 * @return Inverse logit applied to each value in x.
 */
template <
    typename T, require_not_var_matrix_t<T>* = nullptr,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr>
inline auto inv_logit(const T& x) {
  return apply_scalar_unary<inv_logit_fun, T>::apply(x);
}


} // namespace math
} // namespace stan
#endif 

