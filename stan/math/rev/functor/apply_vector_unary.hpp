#ifndef STAN_MATH_REV_FUNCTOR_APPLY_VECTOR_UNARY_HPP
#define STAN_MATH_REV_FUNCTOR_APPLY_VECTOR_UNARY_HPP

#include <stan/math/prim/functor/apply_vector_unary.hpp>
#include <stan/math/rev/core/var.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Specialisation for use with `var_value<T>` types where T inherits from
 * EigenBase. Inputs are passed through unmodified.
 */
template <typename T>
struct apply_vector_unary<T, require_var_matrix_t<T>> {
  /**
   * Member function for applying a functor to a `var_value<T>` and
   * subsequently returning a `var_value<T>`.
   *
   * @tparam T Type of argument to which functor is applied.
   * @tparam F Type of functor to apply.
   * @param x input to which operation is applied.
   * @param f functor to apply to Eigen input.
   * @return object with result of applying functor to input
   */
  template <typename F>
  static inline plain_type_t<T> apply(const T& x, const F& f) {
    return f(x);
  }

  /**
   * Member function for applying a functor to a `var_value<T>` and
   * subsequently returning a `var_value<T>`.
   *
   * @tparam T Type of argument to which functor is applied.
   * @tparam F Type of functor to apply.
   * @param x input to which operation is applied.
   * @param f functor to apply to Eigen input.
   * @return object with result of applying functor to input
   */
  template <typename F>
  static inline plain_type_t<T> apply_no_holder(const T& x, const F& f) {
    return f(x);
  }

  /**
   * Member function for applying a functor to a `var_value<T>` and
   * subsequently returning a var. The reduction to a var needs
   * to be implemented in the definition of the functor.
   *
   * @tparam T Type of argument to which functor is applied.
   * @tparam F Type of functor to apply.
   * @param x input to which operation is applied.
   * @param f functor to apply to input.
   * @return scalar result of applying functor to input.
   */
  template <typename F>
  static inline var reduce(const T& x, const F& f) {
    return f(x);
  }
};

}  // namespace math
}  // namespace stan
#endif
