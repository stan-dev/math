#ifndef STAN_MATH_OPENCL_REV_MULTIPLY_LOG_HPP
#define STAN_MATH_OPENCL_REV_MULTIPLY_LOG_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/rev/adjoint_results.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/adjoint_of.hpp>
#include <stan/math/rev/fun/value_of.hpp>

namespace stan {
namespace math {
namespace internalcl {
template <bool Cond, typename T>
inline decltype(auto) conditional_sum(T&& x) {
  if constexpr (Cond) {
    return sum(std::forward<T>(x));
  } else {
    return std::forward<T>(x);
  }
}
}  // namespace internalcl

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
          require_all_prim_or_rev_kernel_expression_t<T_a, T_b>* = nullptr,
          require_any_var_t<T_a, T_b>* = nullptr,
          require_any_not_stan_scalar_t<T_a, T_b>* = nullptr>
inline var_value<matrix_cl<double>> multiply_log(T_a&& a, T_b&& b) {
  arena_t<T_a> a_arena = std::forward<T_a>(a);
  arena_t<T_b> b_arena = std::forward<T_b>(b);

  return make_callback_var(
      multiply_log(value_of(a_arena), value_of(b_arena)),
      [a_arena, b_arena](const vari_value<matrix_cl<double>>& res) mutable {
        constexpr bool is_scalar_a = !is_matrix_v<T_a>;
        constexpr bool is_scalar_b = !is_matrix_v<T_b>;
        using internalcl::conditional_sum;
        auto is_zero = value_of(a_arena) == 0.0 && value_of(b_arena) == 0.0;
        if constexpr (is_var<T_a>::value && is_var<T_b>::value) {
          a_arena.adj() += conditional_sum<is_scalar_a>(select(
              is_zero, 0.0, elt_multiply(res.adj(), log(value_of(b_arena)))));
          b_arena.adj() += conditional_sum<is_scalar_b>(select(
              is_zero, 0.0,
              elt_multiply(res.adj(),
                           elt_divide(value_of(a_arena), value_of(b_arena)))));
        } else if constexpr (is_var<T_a>::value) {
          a_arena.adj() += conditional_sum<is_scalar_a>(select(
              is_zero, 0.0, elt_multiply(res.adj(), log(value_of(b_arena)))));
        } else if constexpr (is_var<T_b>::value) {
          b_arena.adj() += conditional_sum<is_scalar_b>(select(
              is_zero, 0.0,
              elt_multiply(res.adj(),
                           elt_divide(value_of(a_arena), value_of(b_arena)))));
        }
      });
}

}  // namespace math
}  // namespace stan

#endif
#endif
