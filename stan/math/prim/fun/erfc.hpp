#ifndef STAN_MATH_PRIM_FUN_ERFC_HPP
#define STAN_MATH_PRIM_FUN_ERFC_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Structure to wrap the `erfc()`
 * so that it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return Complementary error function applied to x.
 */
struct erfc_fun {
  template <typename T>
  static inline T fun(const T& x) {
    using std::erfc;
    return erfc(x);
  }
};

/**
 * Returns the elementwise `erfc()` of the input,
 * which may be a scalar or any Stan container of numeric scalars.
 *
 * @tparam T type of container
 * @param x container
 * @return Complementary error function applied to each value in x.
 */
template <
    typename T,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr,
    require_not_var_matrix_t<T>* = nullptr>
inline auto erfc(const T& x) {
  return apply_scalar_unary<erfc_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
