
#ifndef STAN_MATH_VECTORIZED_FUN_TGAMMA_HPP
#define STAN_MATH_VECTORIZED_FUN_TGAMMA_HPP
#include <stan/math/prim/fun/tgamma.hpp>
namespace stan {
namespace math {

/**
 * Vectorized version of tgamma().
 *
 * @tparam T type of container
 * @param x container
 * @return Gamma function applied to each value in x.
 * @throw std::domain_error if any value is 0 or a negative integer
 */
template <
    typename T,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr,
    require_not_var_matrix_t<T>* = nullptr>
inline auto tgamma(const T& x) {
  return apply_scalar_unary<tgamma_fun, T>::apply(x);
}


} // namespace math
} // namespace stan
#endif 

