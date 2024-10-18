
#ifndef STAN_MATH_VECTORIZED_FUN_EXP2_HPP
#define STAN_MATH_VECTORIZED_FUN_EXP2_HPP
#include <stan/math/prim/fun/exp2.hpp>
namespace stan {
namespace math {

/**
 * Return the elementwise `exp2()` of the specified argument,
 * which may be a scalar or any Stan container of numeric scalars.
 * The return type is the same as the argument type.
 *
 * @tparam T type of container
 * @param x container
 * @return Elementwise exp2 of members of container.
 */
template <
    typename T,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr,
    require_not_var_matrix_t<T>* = nullptr>
inline auto exp2(const T& x) {
  return apply_scalar_unary<exp2_fun, T>::apply(x);
}


} // namespace math
} // namespace stan
#endif 

