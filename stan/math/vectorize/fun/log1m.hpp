
#ifndef STAN_MATH_VECTORIZED_FUN_LOG1M_HPP
#define STAN_MATH_VECTORIZED_FUN_LOG1M_HPP
#include <stan/math/prim/fun/log1m.hpp>
namespace stan {
namespace math {

/**
 * Vectorized version of log1m().
 *
 * @tparam T type of container
 * @param x container
 * @return Natural log of 1 minus each value in x.
 */
template <
    typename T, require_not_var_matrix_t<T>* = nullptr,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr>
inline auto log1m(const T& x) {
  return apply_scalar_unary<log1m_fun, T>::apply(x);
}


} // namespace math
} // namespace stan
#endif 

