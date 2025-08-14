#ifndef STAN_MATH_PRIM_FUN_EXP2_HPP
#define STAN_MATH_PRIM_FUN_EXP2_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T, require_arithmetic_t<T>* = nullptr>
inline auto exp2(T&& x) {
  return std::exp2(x);
}
/**
 * Structure to wrap `exp2()` so it can be vectorized.
 */
struct exp2_fun {
  /**
   * Return the base two exponent of the specified argument.
   *
   * @tparam T type of argument
   * @param x argument
   * @return Base two exponent of the argument.
   */
  template <typename T>
  static inline auto fun(T&& x) {
    return exp2(std::forward<T>(x));
  }
};

/**
 * Return the elementwise `exp2()` of the specified argument,
 * which may be a scalar or any Stan container of numeric scalars.
 * The return type is the same as the argument type.
 *
 * @tparam T type of container
 * @param x container
 * @return Elementwise exp2 of members of container.
 */
template <
    typename T,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr,
    require_container_t<T>* = nullptr, require_not_var_matrix_t<T>* = nullptr>
inline auto exp2(T&& x) {
  return apply_scalar_unary<exp2_fun, T>::apply(std::forward<T>(x));
}

}  // namespace math
}  // namespace stan

#endif
