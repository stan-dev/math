#ifndef STAN_MATH_PRIM_FUN_CBRT_HPP
#define STAN_MATH_PRIM_FUN_CBRT_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Structure to wrap `cbrt()` so it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return Cube root of x.
 */
struct cbrt_fun {
  template <typename T>
  static inline T fun(const T& x) {
    using std::cbrt;
    return cbrt(x);
  }
};

/**
 * Returns the elementwise `cbrt()` of the input,
 * which may be a scalar or any Stan container of numeric scalars.
 *
 * @tparam T type of container
 * @param x container
 * @return Cube root of each value in x.
 */
template <
    typename T, require_not_var_matrix_t<T>* = nullptr,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr>
inline auto cbrt(const T& x) {
  return apply_scalar_unary<cbrt_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
