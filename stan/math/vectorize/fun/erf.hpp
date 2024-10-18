
#ifndef STAN_MATH_VECTORIZED_FUN_ERF_HPP
#define STAN_MATH_VECTORIZED_FUN_ERF_HPP
#include <stan/math/prim/fun/erf.hpp>
namespace stan {
namespace math {

/**
 * Returns the elementwise `erf()` of the input,
 * which may be a scalar or any Stan container of numeric scalars.
 *
 * @tparam T type of container
 * @param x container
 * @return Error function applied to each value in x.
 */
template <
    typename T,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr,
    require_not_var_matrix_t<T>* = nullptr>
inline auto erf(const T& x) {
  return apply_scalar_unary<erf_fun, T>::apply(x);
}


} // namespace math
} // namespace stan
#endif 

