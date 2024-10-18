
#ifndef STAN_MATH_VECTORIZED_FUN_CBRT_HPP
#define STAN_MATH_VECTORIZED_FUN_CBRT_HPP
#include <stan/math/prim/fun/cbrt.hpp>
namespace stan {
namespace math {

/**
 * Returns the elementwise `cbrt()` of the input,
 * which may be a scalar or any Stan container of numeric scalars.
 *
 * @tparam T type of container
 * @param x container
 * @return Cube root of each value in x.
 */
template <
    typename T, require_not_var_matrix_t<T>* = nullptr,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr>
inline auto cbrt(const T& x) {
  return apply_scalar_unary<cbrt_fun, T>::apply(x);
}


} // namespace math
} // namespace stan
#endif 

