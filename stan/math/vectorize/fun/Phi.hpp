
#ifndef STAN_MATH_VECTORIZED_FUN_PHI_HPP
#define STAN_MATH_VECTORIZED_FUN_PHI_HPP
#include <stan/math/prim/fun/Phi.hpp>
namespace stan {
namespace math {

/**
 * Vectorized version of Phi().
 *
 * @tparam T type of container
 * @param x container
 * @return Unit normal CDF of each value in x.
 */
template <
    typename T,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr,
    require_not_var_matrix_t<T>* = nullptr>
inline auto Phi(const T& x) {
  return apply_scalar_unary<Phi_fun, T>::apply(x);
}


} // namespace math
} // namespace stan
#endif 

