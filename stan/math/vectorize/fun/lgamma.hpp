
#ifndef STAN_MATH_VECTORIZED_FUN_LGAMMA_HPP
#define STAN_MATH_VECTORIZED_FUN_LGAMMA_HPP
#include <stan/math/prim/fun/lgamma.hpp>
namespace stan {
namespace math {

/**
 * Vectorized version of lgamma().
 *
 * @tparam T type of container
 * @param x container
 * @return Natural log of the gamma function
 *         applied to each value in x.
 * @throw std::domain_error if any value is a negative integer or 0.
 */
template <typename T, require_not_var_matrix_t<T>* = nullptr,
          require_not_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr>
inline auto lgamma(const T& x) {
  return apply_scalar_unary<lgamma_fun, T>::apply(x);
}


} // namespace math
} // namespace stan
#endif 

