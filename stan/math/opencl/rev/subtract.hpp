#ifndef STAN_MATH_OPENCL_REV_SUBTRACT_HPP
#define STAN_MATH_OPENCL_REV_SUBTRACT_HPP
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
inline auto subtract(const T_a& a, const T_b& b) {
  check_size_match("subtract (OpenCL)", "A.cols()", a.cols(), "B.cols()",
                   b.cols());
  check_size_match("subtract (OpenCL)", "A.rows()", a.rows(), "B.rows()",
                   b.rows());
  const arena_t<T_a>& a_arena = a;
  const arena_t<T_b>& b_arena = b;

  var_value<matrix_cl<double>> res = value_of(a_arena) - value_of(b_arena);

  reverse_pass_callback([a_arena, b_arena, res]() mutable {
    if (!is_constant<T_a>::value && !is_constant<T_b>::value) {
      auto& a_adj = forward_as<var_value<matrix_cl<double>>>(a_arena).adj();
      auto& b_adj = forward_as<var_value<matrix_cl<double>>>(b_arena).adj();
      results(a_adj, b_adj)
          = expressions((a_adj + res.adj()), (b_adj - res.adj()));
    } else if (!is_constant<T_a>::value) {
      auto& a_adj = forward_as<var_value<matrix_cl<double>>>(a_arena).adj();
      a_adj = a_adj + res.adj();
    } else {
      auto& b_adj = forward_as<var_value<matrix_cl<double>>>(b_arena).adj();
      b_adj = b_adj - res.adj();
    }
  });
  return res;
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

/**
 * Subtraction of a kernel generator expression and a scalar.
 *
 * @tparam T1 type of the scalar
 * @tparam T2 type of the matrix or expression
 *
 * @param a scalar
 * @param b matrix
 * @return The subtraction of the the second argument from the first
 */
template <typename T1, typename T2, require_stan_scalar_t<T1>* = nullptr,
          require_all_nonscalar_prim_or_rev_kernel_expression_t<T2>* = nullptr,
          require_any_var_t<T1, T2>* = nullptr>
inline auto subtract(const T1& a, const T2& b) {
  const arena_t<T1>& a_arena = a;
  const arena_t<T2>& b_arena = b;

  var_value<matrix_cl<double>> res = value_of(a_arena) - value_of(b_arena);

  reverse_pass_callback([a_arena, b_arena, res]() mutable {
    if (!is_constant<T1>::value) {
      auto& a_adj = forward_as<var_value<double>>(a_arena).adj();
      a_adj = a_adj + sum(res.adj());
    }
    if (!is_constant<T2>::value) {
      auto& b_adj = forward_as<var_value<matrix_cl<double>>>(b_arena).adj();
      b_adj = b_adj - res.adj();
    }
  });
  return res;
}

/**
 * Subtraction of a kernel generator expression and a scalar.
 *
 * @tparam T1 type of the matrix or expression
 * @tparam T2 type of the scalar
 *
 * @param a matrix
 * @param b scalar
 * @return The subtraction of the the second argument from the first
 */
template <typename T1, typename T2, require_stan_scalar_t<T2>* = nullptr,
          require_all_nonscalar_prim_or_rev_kernel_expression_t<T1>* = nullptr,
          require_any_var_t<T1, T2>* = nullptr>
inline auto subtract(const T1& a, const T2& b) {
  const arena_t<T1>& a_arena = a;
  const arena_t<T2>& b_arena = b;

  var_value<matrix_cl<double>> res = value_of(a_arena) - value_of(b_arena);

  reverse_pass_callback([a_arena, b_arena, res]() mutable {
    if (!is_constant<T1>::value) {
      auto& a_adj = forward_as<var_value<matrix_cl<double>>>(a_arena).adj();
      a_adj = a_adj + res.adj();
    }
    if (!is_constant<T2>::value) {
      auto& b_adj = forward_as<var_value<double>>(b_arena).adj();
      b_adj = b_adj - sum(res.adj());
    }
  });
  return res;
}

}  // namespace math
}  // namespace stan

#endif
#endif
