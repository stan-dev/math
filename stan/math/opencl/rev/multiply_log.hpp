#ifndef STAN_MATH_OPENCL_REV_MULTIPLY_LOG_HPP
#define STAN_MATH_OPENCL_REV_MULTIPLY_LOG_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/adjoint_of.hpp>
#include <stan/math/rev/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Returns the elementwise `multiply_log()` of the input.
 *
 * @tparam T_a type of first expression
 * @tparam T_b type of second expression
 * @param a first expression
 * @param b second expression
 *
 * @return Elementwise `multiply_log()` of the input.
 */
template <
    typename T_a, typename T_b,
    require_all_nonscalar_prim_or_rev_kernel_expression_t<T_a, T_b>* = nullptr,
    require_any_var_t<T_a, T_b>* = nullptr>
inline var_value<matrix_cl<double>> multiply_log(T_a&& a, T_b&& b) {
  const arena_t<T_a>& a_arena = std::forward<T_a>(a);
  const arena_t<T_b>& b_arena = std::forward<T_b>(b);

  matrix_cl<double> res_val
      = multiply_log(value_of(a_arena), value_of(b_arena));

  return make_callback_var(
      res_val,
      [a_arena, b_arena](const vari_value<matrix_cl<double>>& res) mutable {
        auto a_deriv = elt_multiply(res.adj(), log(value_of(b_arena)));
        auto b_deriv = elt_multiply(
            res.adj(), elt_divide(value_of(a_arena), value_of(b_arena)));
        results(adjoint_of(a_arena), adjoint_of(b_arena))
            += expressions(calc_if<is_var<T_a>::value>(a_deriv),
                           calc_if<is_var<T_b>::value>(b_deriv));
      });
}

/**
 * Returns the elementwise `multiply_log()` of the input.
 *
 * @tparam T_a type of first expression
 * @tparam T_b type of second expression
 * @param a first expression
 * @param b second expression
 *
 * @return Elementwise `multiply_log()` of the input.
 */
template <typename T_a, typename T_b,
          require_nonscalar_prim_or_rev_kernel_expression_t<T_a>* = nullptr,
          require_stan_scalar_t<T_b>* = nullptr,
          require_any_var_t<T_a, T_b>* = nullptr>
inline var_value<matrix_cl<double>> multiply_log(T_a&& a, const T_b& b) {
  const arena_t<T_a>& a_arena = std::forward<T_a>(a);

  matrix_cl<double> res_val = multiply_log(value_of(a_arena), value_of(b));

  return make_callback_var(
      res_val, [a_arena, b](const vari_value<matrix_cl<double>>& res) mutable {
        if (!is_constant<T_a>::value) {
          auto& a_adj = forward_as<var_value<matrix_cl<double>>>(a_arena).adj();
          a_adj = a_adj + elt_multiply(res.adj(), log(value_of(b)));
        }
        if (!is_constant<T_b>::value) {
          adjoint_of(b) += sum(elt_multiply(
              res.adj(), elt_divide(value_of(a_arena), value_of(b))));
        }
      });
}
/**
 * Returns the elementwise `multiply_log()` of the input.
 *
 * @tparam T_a type of first expression
 * @tparam T_b type of second expression
 * @param a first expression
 * @param b second expression
 *
 * @return Elementwise `multiply_log()` of the input.
 */
template <typename T_a, typename T_b,
          require_nonscalar_prim_or_rev_kernel_expression_t<T_b>* = nullptr,
          require_stan_scalar_t<T_a>* = nullptr,
          require_any_var_t<T_a, T_b>* = nullptr>
inline var_value<matrix_cl<double>> multiply_log(const T_a& a, T_b&& b) {
  const arena_t<T_b>& b_arena = std::forward<T_b>(b);

  matrix_cl<double> res_val = multiply_log(value_of(a), value_of(b_arena));

  return make_callback_var(
      res_val, [a, b_arena](const vari_value<matrix_cl<double>>& res) mutable {
        if (!is_constant<T_a>::value) {
          adjoint_of(a) += sum(elt_multiply(res.adj(), log(value_of(b_arena))));
        }
        if (!is_constant<T_b>::value) {
          auto& b_adj = forward_as<var_value<matrix_cl<double>>>(b_arena).adj();
          b_adj = b_adj
                  + elt_multiply(res.adj(),
                                 elt_divide(value_of(a), value_of(b_arena)));
        }
      });
}

}  // namespace math
}  // namespace stan

#endif
#endif
