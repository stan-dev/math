#ifndef STAN_MATH_OPENCL_REV_HYPOT_HPP
#define STAN_MATH_OPENCL_REV_HYPOT_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/adjoint_of.hpp>
#include <stan/math/rev/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Returns the elementwise `hypot()` of the input.
 *
 * @tparam T_a type of first expression
 * @tparam T_b type of second expression
 * @param a first expression
 * @param b second expression
 *
 * @return Elementwise `hypot()` of the input.
 */
template <
    typename T_a, typename T_b,
    require_all_nonscalar_prim_or_rev_kernel_expression_t<T_a, T_b>* = nullptr,
    require_any_var_t<T_a, T_b>* = nullptr>
inline var_value<matrix_cl<double>> hypot(T_a&& a, T_b&& b) {
  const arena_t<T_a>& a_arena = std::forward<T_a>(a);
  const arena_t<T_b>& b_arena = std::forward<T_b>(b);

  matrix_cl<double> res_val = hypot(value_of(a_arena), value_of(b_arena));

  return make_callback_var(
      res_val,
      [a_arena, b_arena](const vari_value<matrix_cl<double>>& res) mutable {
        auto a_deriv
            = elt_multiply(res.adj(), elt_divide(value_of(a_arena), res.val()));
        auto b_deriv
            = elt_multiply(res.adj(), elt_divide(value_of(b_arena), res.val()));
        results(adjoint_of(a_arena), adjoint_of(b_arena))
            += expressions(calc_if<is_var<T_a>::value>(a_deriv),
                           calc_if<is_var<T_b>::value>(b_deriv));
      });
}

/**
 * Returns the elementwise `hypot()` of the input, where the first input
 * is a nonscalar kernel generator expression and the second input is
 * a scalar.
 *
 * @tparam T_a type of kernel generator expression
 * @tparam T_b type of scalar
 * @param a kernel generator expression
 * @param b scalar
 * @return Elementwise `hypot()` of the input arguments
 */
template <typename T_a, typename T_b,
          require_nonscalar_prim_or_rev_kernel_expression_t<T_a>* = nullptr,
          require_stan_scalar_t<T_b>* = nullptr,
          require_any_var_t<T_a, T_b>* = nullptr>
inline var_value<matrix_cl<double>> hypot(T_a&& a, const T_b& b) {
  const arena_t<T_a>& a_arena = std::forward<T_a>(a);

  matrix_cl<double> res_val = hypot(value_of(a_arena), value_of(b));

  return make_callback_var(
      res_val, [a_arena, b](const vari_value<matrix_cl<double>>& res) mutable {
        if (!is_constant<T_a>::value) {
          auto& a_adj = forward_as<var_value<matrix_cl<double>>>(a_arena).adj();
          a_adj = a_adj
                  + elt_multiply(res.adj(),
                                 elt_divide(value_of(a_arena), res.val()));
        }
        if (!is_constant<T_b>::value) {
          adjoint_of(b) += sum(
              elt_multiply(res.adj(), elt_divide(value_of(b), res.val())));
        }
      });
}

/**
 * Returns the elementwise `hypot()` of the input, where the first input
 * is a scalar and the second input is a nonscalar kernel generator
 * expression.
 *
 * @tparam T_a type of scalar
 * @tparam T_b type of kernel generator expression
 * @param a scalar
 * @param b kernel generator expression
 * @return Elementwise division of the input arguments
 */
template <typename T_a, typename T_b,
          require_nonscalar_prim_or_rev_kernel_expression_t<T_b>* = nullptr,
          require_stan_scalar_t<T_a>* = nullptr,
          require_any_var_t<T_a, T_b>* = nullptr>
inline var_value<matrix_cl<double>> hypot(const T_a& a, T_b&& b) {
  const arena_t<T_b>& b_arena = std::forward<T_b>(b);

  matrix_cl<double> res_val = hypot(value_of(a), value_of(b_arena));

  return make_callback_var(
      res_val, [a, b_arena](const vari_value<matrix_cl<double>>& res) mutable {
        if (!is_constant<T_a>::value) {
          adjoint_of(a) += sum(
              elt_multiply(res.adj(), elt_divide(value_of(a), res.val())));
        }
        if (!is_constant<T_b>::value) {
          auto& b_adj = forward_as<var_value<matrix_cl<double>>>(b_arena).adj();
          b_adj = b_adj
                  + elt_multiply(res.adj(),
                                 elt_divide(value_of(b_arena), res.val()));
        }
      });
}

}  // namespace math
}  // namespace stan

#endif
#endif
