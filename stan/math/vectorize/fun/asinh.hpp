
#ifndef STAN_MATH_VECTORIZED_FUN_ASINH_HPP
#define STAN_MATH_VECTORIZED_FUN_ASINH_HPP
#include <stan/math/prim/fun/asinh.hpp>
namespace stan {
namespace math {

/**
 * Returns the elementwise `asinh()` of the input,
 * which may be a scalar or any Stan container of numeric scalars.
 *
 * @tparam T type of container
 * @param x container
 * @return Inverse hyperbolic sine of each value in the container.
 */
template <
    typename T, require_not_var_matrix_t<T>* = nullptr,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr>
inline auto asinh(const T& x) {
  return apply_scalar_unary<asinh_fun, T>::apply(x);
}


} // namespace math
} // namespace stan
#endif 

