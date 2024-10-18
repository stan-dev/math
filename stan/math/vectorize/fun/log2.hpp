
#ifndef STAN_MATH_VECTORIZED_FUN_LOG2_HPP
#define STAN_MATH_VECTORIZED_FUN_LOG2_HPP
#include <stan/math/prim/fun/log2.hpp>
namespace stan {
namespace math {

/**
 * Return the elementwise application of `log2()` to
 * specified argument container.  The return type promotes the
 * underlying scalar argument type to double if it is an integer,
 * and otherwise is the argument type.
 *
 * @tparam T type of container
 * @param x container
 * @return elementwise log2 of container elements
 */
template <typename T, require_not_var_matrix_t<T>* = nullptr,
          require_not_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr>
inline auto log2(const T& x) {
  return apply_scalar_unary<log2_fun, T>::apply(x);
}


} // namespace math
} // namespace stan
#endif 

