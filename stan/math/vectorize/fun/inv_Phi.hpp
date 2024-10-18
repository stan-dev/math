
#ifndef STAN_MATH_VECTORIZED_FUN_INV_PHI_HPP
#define STAN_MATH_VECTORIZED_FUN_INV_PHI_HPP
#include <stan/math/prim/fun/inv_Phi.hpp>
namespace stan {
namespace math {

/**
 * Vectorized version of inv_Phi().
 *
 * @tparam T type of container
 * @param x variables in range [0, 1]
 * @return Inverse unit normal CDF of each value in x.
 * @throw std::domain_error if any value is not between 0 and 1.
 */
template <
    typename T,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr,
    require_not_var_matrix_t<T>* = nullptr>
inline auto inv_Phi(const T& x) {
  return apply_scalar_unary<inv_Phi_fun, T>::apply(x);
}


} // namespace math
} // namespace stan
#endif 

