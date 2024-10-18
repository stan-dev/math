
#ifndef STAN_MATH_VECTORIZED_FUN_PHI_APPROX_HPP
#define STAN_MATH_VECTORIZED_FUN_PHI_APPROX_HPP
#include <stan/math/prim/fun/Phi_approx.hpp>
namespace stan {
namespace math {

/**
 * Return the elementwise application of <code>Phi_approx()</code> to
 * specified argument container.  The return type promotes the
 * underlying scalar argument type to double if it is an integer,
 * and otherwise is the argument type.
 *
 * @tparam T type of container
 * @param x container
 * @return elementwise Phi_approx of container elements
 */
template <
    typename T,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr,
    require_not_var_matrix_t<T>* = nullptr>
inline auto Phi_approx(const T& x) {
  return apply_scalar_unary<Phi_approx_fun, T>::apply(x);
}


} // namespace math
} // namespace stan
#endif 

