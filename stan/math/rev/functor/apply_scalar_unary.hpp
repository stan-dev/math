#ifndef STAN_MATH_REV_FUNCTOR_APPLY_SCALAR_UNARY_HPP
#define STAN_MATH_REV_FUNCTOR_APPLY_SCALAR_UNARY_HPP

#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <stan/math/rev/core/var.hpp>

namespace stan {
namespace math {

/**
 * Template specialization to var for vectorizing a unary scalar
 * function.  This is a base scalar specialization.  It applies
 * the function specified by the template parameter to the
 * argument.
 *
 * @tparam F Type of function to apply.
 */
template <typename F, typename T>
struct apply_scalar_unary<
    F, T, std::enable_if_t<is_var<T>::value && is_stan_scalar<T>::value>> {
  /**
   * Function return type, which is <code>var</code>.
   */
  using return_t = var;

  /**
   * Apply the function specified by F to the specified argument.
   *
   * @param x Argument variable.
   * @return Function applied to the variable.
   */
  template <typename T2>
  static inline auto apply(T2&& x) {
    return F::fun(std::forward<T2>(x));
  }
};

template <typename F, typename T>
struct apply_scalar_unary<F, T, require_var_matrix_t<T>> {
  /**
   * Function return type, which is a `var_value` with plain value type.
   */
  using return_t = plain_type_t<T>;

  /**
   * Apply the function specified by F to the specified argument.
   *
   * @param x Argument variable.
   * @return Function applied to the variable.
   */
  template <typename T2>
  static inline auto apply(T2&& x) {
    return F::fun(std::forward<T2>(x));
  }
};

}  // namespace math
}  // namespace stan
#endif
