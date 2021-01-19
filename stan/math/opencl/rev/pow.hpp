#ifndef STAN_MATH_OPENCL_REV_POW_HPP
#define STAN_MATH_OPENCL_REV_POW_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/rev/adjoint_results.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/sum.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/adjoint_of.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/rev/core/reverse_pass_callback.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/meta/is_kernel_expression.hpp>

namespace stan {
namespace math {

/**
 * Return the first argument raised to the power of the second
 * argument.
 *
 * @tparam T_a type of first expression of the base
 * @tparam T_b type of second expression of the exponent
 * @param a first expression of base
 * @param b second expression of exponent
 * @return base raised to the power of the exponent
 */
template <typename T_a, typename T_b,
          require_all_prim_or_rev_kernel_expression_t<T_a, T_b>* = nullptr,
          require_any_var_t<T_a, T_b>* = nullptr,
          require_any_not_stan_scalar_t<T_a, T_b>* = nullptr>
inline var_value<matrix_cl<double>> pow(T_a&& a, T_b&& b) {
  const arena_t<T_a>& a_arena = std::forward<T_a>(a);
  const arena_t<T_b>& b_arena = std::forward<T_b>(b);

  matrix_cl<double> res_val = pow(value_of(a_arena), value_of(b_arena));

  return make_callback_var(
      res_val,
      [a_arena, b_arena](const vari_value<matrix_cl<double>>& res) mutable {
        auto isnan_res = isnan(res.val());
        auto zero_a_arena = value_of(a_arena) == 0.0;
        auto res_adj_times_val = elt_multiply(res.adj(), res.val());
        auto a_deriv = select(
            isnan_res, NOT_A_NUMBER,
            select(zero_a_arena, 0.0,
                   elt_multiply(
                       res_adj_times_val,
                       elt_divide(value_of(b_arena), value_of(a_arena)))));
        auto b_deriv = select(
            isnan_res, NOT_A_NUMBER,
            select(zero_a_arena, 0.0,
                   elt_multiply(res_adj_times_val, log(value_of(a_arena)))));

        adjoint_results(a_arena, b_arena) += expressions(a_deriv, b_deriv);
      });
}

}  // namespace math
}  // namespace stan

#endif
#endif
