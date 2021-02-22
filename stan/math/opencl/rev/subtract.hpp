#ifndef STAN_MATH_OPENCL_REV_SUBTRACT_HPP
#define STAN_MATH_OPENCL_REV_SUBTRACT_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/rev/adjoint_results.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/sum.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/rev/core/reverse_pass_callback.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/meta/is_kernel_expression.hpp>

namespace stan {
namespace math {

/**
 * Subtraction of two reverse mode matrices and/or kernel generator
 * expressions.
 * @tparam T_a type of first expression
 * @tparam T_b type of second expression
 * @param a first expression
 * @param b second expression
 * @return The subtraction of the the second argument from the first
 */
template <typename T_a, typename T_b,
          require_all_prim_or_rev_kernel_expression_t<T_a, T_b>* = nullptr,
          require_any_var_t<T_a, T_b>* = nullptr,
          require_any_not_stan_scalar_t<T_a, T_b>* = nullptr>
inline auto subtract(T_a&& a, T_b&& b) {
  arena_t<T_a> a_arena = a;
  arena_t<T_b> b_arena = b;

  return make_callback_var(
      value_of(a_arena) - value_of(b_arena),
      [a_arena, b_arena](vari_value<matrix_cl<double>>& res) mutable {
        adjoint_results(a_arena, b_arena) += expressions(res.adj(), -res.adj());
      });
}

/**
 * Subtraction of two reverse mode matrices and/or kernel generator
 * expressions.
 * @tparam T_a type of first expression
 * @tparam T_b type of second expression
 * @param a first expression
 * @param b second expression
 * @return The subtraction of the the second argument from the first
 */
template <
    typename T_a, typename T_b,
    require_all_nonscalar_prim_or_rev_kernel_expression_t<T_a, T_b>* = nullptr,
    require_any_var_t<T_a, T_b>* = nullptr>
inline auto operator-(const T_a& a, const T_b& b) {
  return subtract(a, b);
}

}  // namespace math
}  // namespace stan

#endif
#endif
