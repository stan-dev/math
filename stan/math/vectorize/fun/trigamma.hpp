
#ifndef STAN_MATH_VECTORIZED_FUN_TRIGAMMA_HPP
#define STAN_MATH_VECTORIZED_FUN_TRIGAMMA_HPP
#include <stan/math/prim/fun/trigamma.hpp>
namespace stan {
namespace math {

/**
 * Return the elementwise application of `trigamma()` to
 * specified argument container.  The return type promotes the
 * underlying scalar argument type to double if it is an integer,
 * and otherwise is the argument type.
 *
 * @tparam T type of container
 * @param x container
 * @return elementwise trigamma of container elements
 */
template <typename T,
          require_not_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr>
inline auto trigamma(const T& x) {
  return apply_scalar_unary<trigamma_fun, T>::apply(x);
}


} // namespace math
} // namespace stan
#endif 

