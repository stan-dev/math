#ifndef STAN_MATH_OPENCL_REV_APPEND_ROW_HPP
#define STAN_MATH_OPENCL_REV_APPEND_ROW_HPP
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
 * Return the result of stacking the rows of the first argument
 * matrix on top of the second argument matrix.
 *
 * Given input types result in following outputs:
 * (matrix, matrix) -> matrix,
 * (matrix, row_vector) -> matrix,
 * (row_vector, matrix) -> matrix,
 * (row_vector, row_vector) -> matrix,
 * (vector, vector) -> vector.
 *
 * @tparam T_a type of the first matrix
 * @tparam T_b type of the second matrix
 *
 * @param a First matrix.
 * @param b Second matrix.
 * @return Result of stacking first matrix on top of second.
 */
template <typename T_a, typename T_b,
          require_all_prim_or_rev_kernel_expression_t<T_a, T_b>* = nullptr,
          require_any_var_t<T_a, T_b>* = nullptr,
          require_any_not_stan_scalar_t<T_a, T_b>* = nullptr>
inline var_value<matrix_cl<double>> append_row(T_a&& a, T_b&& b) {
  arena_t<T_a> a_arena = std::forward<T_a>(a);
  arena_t<T_b> b_arena = std::forward<T_b>(b);

  return make_callback_var(
      append_row(value_of(a_arena), value_of(b_arena)),
      [a_arena, b_arena](const vari_value<matrix_cl<double>>& res) mutable {
        if (!is_constant<T_a>::value) {
          adjoint_of(a_arena) += block_zero_based(
              res.adj(), 0, 0, a_arena.rows(), a_arena.cols());
        }
        if (!is_constant<T_b>::value) {
          adjoint_of(b_arena) += block_zero_based(
              res.adj(), a_arena.rows(), 0, b_arena.rows(), b_arena.cols());
        }
      });
}

}  // namespace math
}  // namespace stan

#endif
#endif
