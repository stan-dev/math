
#ifndef STAN_MATH_VECTORIZED_FUN_LAMBERT_WM1_HPP
#define STAN_MATH_VECTORIZED_FUN_LAMBERT_WM1_HPP
#include <stan/math/prim/fun/lambert_w.hpp>
namespace stan {
namespace math {

/**
 * Vectorized version of lambert_wm1().
 *
 * @tparam T type of container
 * @param x container
 * @return value of the W0 branch of the Lambert W function for each value in x
 * @throw std::domain_error if x is less than or equal to `-e^(-1)` or greater
 * than or equal to 0
 */
template <typename T, require_not_stan_scalar_t<T>* = nullptr,
          require_not_var_matrix_t<T>* = nullptr>
inline auto lambert_wm1(const T& x) {
  return apply_scalar_unary<internal::lambert_wm1_fun, T>::apply(x);
}


} // namespace math
} // namespace stan
#endif 

