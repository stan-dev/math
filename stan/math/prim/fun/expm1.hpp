#ifndef STAN_MATH_PRIM_FUN_EXPM1_HPP
#define STAN_MATH_PRIM_FUN_EXPM1_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Structure to wrap `expm1()` so that it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return Natural exponential of x minus one.
 */
struct expm1_fun {
  template <typename T>
  static inline T fun(const T& x) {
    using std::expm1;
    return expm1(x);
  }
};

/**
 * Return the elementwise `expm1()` of the specified argument,
 * which may be a scalar or any Stan container of numeric scalars.
 * The return type is the same as the argument type.
 *
 * @tparam T type of container
 * @param x container
 * @return Natural exponential of each value in x minus one.
 */
template <
    typename T,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr,
    require_not_var_matrix_t<T>* = nullptr>
inline auto expm1(const T& x) {
  return apply_scalar_unary<expm1_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
