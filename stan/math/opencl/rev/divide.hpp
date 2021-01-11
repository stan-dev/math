#ifndef STAN_MATH_OPENCL_REV_DIVIDE_HPP
#define STAN_MATH_OPENCL_REV_DIVIDE_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta/is_kernel_expression.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/sum.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/rev/core/reverse_pass_callback.hpp>
#include <stan/math/prim/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Elementwise division of a kernel generator expression and a scalar.
 *
 * @tparam T_a type of the kernel generator expression
 * @tparam T_b type of the scalar
 *
 * @param a input kernel generator expression
 * @param b scalar
 * @return Elementwise division of the input kernel generator
 *  expression with a scalar
 */
template <typename T_a, typename T_b, require_stan_scalar_t<T_b>* = nullptr,
          require_all_nonscalar_prim_or_rev_kernel_expression_t<T_a>* = nullptr,
          require_any_var_t<T_a, T_b>* = nullptr>
inline var_value<matrix_cl<double>> divide(T_a&& a, T_b&& b) {
  arena_t<T_a> a_arena = std::forward<T_a>(a);
  arena_t<T_b> b_arena = std::forward<T_b>(b);

  matrix_cl<double> res_val = elt_divide(value_of(a_arena), value_of(b_arena));

  return make_callback_var(
      res_val,
      [a_arena, b_arena](const vari_value<matrix_cl<double>>& res) mutable {
        if (!is_constant<T_a>::value) {
          auto& a_adj = forward_as<var_value<matrix_cl<double>>>(a_arena).adj();
          a_adj = a_adj + elt_divide(res.adj(), value_of(b_arena));
        }
        if (!is_constant<T_b>::value) {
          auto& b_adj = forward_as<var_value<double>>(b_arena).adj();
          b_adj = b_adj - sum(elt_divide(res.val(), value_of(b_arena)));
        }
      });
}

}  // namespace math
}  // namespace stan

#endif
#endif
