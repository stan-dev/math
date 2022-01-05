#ifndef STAN_MATH_PRIM_FUN_INV_ERFC_HPP
#define STAN_MATH_PRIM_FUN_INV_ERFC_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <boost/math/special_functions/erf.hpp>

namespace stan {
namespace math {

/**
 * The inverse complementary error function for variables.
 *
 * @tparam T type of scalar
 * @param x scalar
 * @return Inverse complementary error function applied to x.
 */
template <typename T, require_arithmetic_t<T>* = nullptr>
inline auto inv_erfc(const T& x) {
  return boost::math::erfc_inv(x);
}

/**
 * Structure to wrap the `inv_erfc() function`
 * so that it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return Inverse of the complementary error function applied to x.
 */
struct inv_erfc_fun {
  template <typename T>
  static inline T fun(const T& x) {
    return inv_erfc(x);
  }
};

/**
 * Returns the elementwise `inv_erfc()` of the input,
 * which may be a scalar or any Stan container of numeric scalars.
 *
 * @tparam T type of container
 * @param x container
 * @return Inverse complementary error function applied to each value in x.
 */
template <
    typename T,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr,
    require_not_var_matrix_t<T>* = nullptr,
    require_not_arithmetic_t<T>* = nullptr>
inline auto inv_erfc(const T& x) {
  return apply_scalar_unary<inv_erfc_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
