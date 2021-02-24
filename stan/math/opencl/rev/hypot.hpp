#ifndef STAN_MATH_OPENCL_REV_HYPOT_HPP
#define STAN_MATH_OPENCL_REV_HYPOT_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/rev/adjoint_results.hpp>
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
template <typename T_a, typename T_b,
          require_all_prim_or_rev_kernel_expression_t<T_a, T_b>* = nullptr,
          require_any_var_t<T_a, T_b>* = nullptr,
          require_any_not_stan_scalar_t<T_a, T_b>* = nullptr>
inline var_value<matrix_cl<double>> hypot(T_a&& a, T_b&& b) {
  arena_t<T_a> a_arena = std::forward<T_a>(a);
  arena_t<T_b> b_arena = std::forward<T_b>(b);

  return make_callback_var(
      hypot(value_of(a_arena), value_of(b_arena)),
      [a_arena, b_arena](const vari_value<matrix_cl<double>>& res) mutable {
        auto res_adj_div_val = elt_divide(res.adj(), res.val());
        adjoint_results(a_arena, b_arena)
            += expressions(elt_multiply(value_of(a_arena), res_adj_div_val),
                           elt_multiply(value_of(b_arena), res_adj_div_val));
      });
}

}  // namespace math
}  // namespace stan

#endif
#endif
