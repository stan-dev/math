
#ifndef STAN_MATH_VECTORIZED_FUN_ATANH_HPP
#define STAN_MATH_VECTORIZED_FUN_ATANH_HPP
#include <stan/math/prim/fun/atanh.hpp>
namespace stan {
namespace math {

/**
 * Return the elementwise application of <code>atanh()</code> to
 * specified argument container.  The return type promotes the
 * underlying scalar argument type to double if it is an integer,
 * and otherwise is the argument type.
 *
 * @tparam T type of container
 * @param x container
 * @return Elementwise atanh of members of container.
 */
template <
    typename T, require_not_var_matrix_t<T>* = nullptr,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr>
inline auto atanh(const T& x) {
  return apply_scalar_unary<atanh_fun, T>::apply(x);
}


} // namespace math
} // namespace stan
#endif 

