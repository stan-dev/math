#ifndef STAN_MATH_OPENCL_REV_MULTIPLY_HPP
#define STAN_MATH_OPENCL_REV_MULTIPLY_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator/is_kernel_expression.hpp>
#include <stan/math/opencl/multiply.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/rev/functor/reverse_pass_callback.hpp>
#include <stan/math/prim/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Matrix multiplication of two reverse mode matrices and/or kernel generator
 * expressions.
 * @tparam T_a type of first expression
 * @tparam T_b type of second expression
 * @param a first expression
 * @param b second expression
 * @return Matrix product of given arguments
 */
template <
    typename T_a, typename T_b,
    require_all_nonscalar_prim_or_rev_kernel_expression_t<T_a, T_b>* = nullptr,
    require_any_var_t<T_a, T_b>* = nullptr>
inline auto multiply(const T_a& a, const T_b& b) {
  check_size_match("multiply ((OpenCL))", "A.cols()", a.cols(), "B.rows()",
                   b.rows());
  const arena_t<T_a>& a_arena = a;
  const arena_t<T_b>& b_arena = b;

  var_value<matrix_cl<double>> res = value_of(a_arena) * value_of(b_arena);

  reverse_pass_callback([a_arena, b_arena, res]() mutable {
    if (!is_constant<T_a>::value) {
      auto& a_adj = forward_as<var_value<matrix_cl<double>>>(a_arena).adj();
      a_adj = a_adj + res.adj() * transpose(value_of(b_arena));
    }
    if (!is_constant<T_b>::value) {
      auto& b_adj = forward_as<var_value<matrix_cl<double>>>(b_arena).adj();
      b_adj = b_adj + transpose(value_of(a_arena)) * res.adj();
    }
  });
  return res;
}

/**
 * Matrix multiplication of two reverse mode matrices and/or kernel generator
 * expressions.
 * @tparam T_a type of first expression
 * @tparam T_b type of second expression
 * @param a first expression
 * @param b second expression
 * @return Matrix product of given arguments
 */
template <
    typename T_a, typename T_b,
    require_all_nonscalar_prim_or_rev_kernel_expression_t<T_a, T_b>* = nullptr,
    require_any_var_t<T_a, T_b>* = nullptr>
inline auto operator*(const T_a& a, const T_b& b) {
  return multiply(a, b);
}

}  // namespace math
}  // namespace stan

#endif
#endif
