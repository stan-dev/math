#ifndef STAN_MATH_OPENCL_REV_MDIVIDE_LEFT_TRI_LOW_HPP
#define STAN_MATH_OPENCL_REV_MDIVIDE_LEFT_TRI_LOW_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/rev/arena_type.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/prim/mdivide_left_tri_low.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/adjoint_of.hpp>
#include <stan/math/rev/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Returns the solution of the system Ax=b when A is lower triangular.
 *
 * @tparam T1 type of elements in A
 * @tparam T2 type of elements in b
 * @param A Triangular matrix.
 * @param b Right hand side matrix or vector.
 * @return x = A^-1 b, solution of the linear system.
 * @throws std::domain_error if A is not square or the rows of b don't
 * match the size of A.
 */
template <
    typename T1, typename T2,
    require_all_nonscalar_prim_or_rev_kernel_expression_t<T1, T2>* = nullptr,
    require_any_var_t<T1, T2>* = nullptr>
inline var_value<matrix_cl<double>> mdivide_left_tri_low(T1&& A, T2&& b) {
  check_square("mdivide_left_tri_low", "A", A);
  check_multiplicable("mdivide_left_tri_low", "A", A, "b", b);
  if (A.size() == 0 || b.size() == 0) {
    return var_value<matrix_cl<double>>(matrix_cl<double>(A.rows(), b.cols()));
  }
  arena_t<T1> A_arena = std::forward<T1>(A);
  arena_t<T2> b_arena = std::forward<T2>(b);
  arena_matrix_cl<double> A_tri_inv
      = tri_inverse<matrix_cl_view::Lower>(value_of(A_arena));
  return make_callback_var(
      A_tri_inv * value_of(b_arena),
      [A_arena, b_arena, A_tri_inv](const vari_value<matrix_cl<double>>& res) {
        matrix_cl<double> adjB = transpose(A_tri_inv) * res.adj();
        if (!is_constant<T1>::value) {
          matrix_cl<double> adjA = adjB * transpose(res.val());
          adjA.view(matrix_cl_view::Lower);
          adjoint_of(A_arena) -= adjA;
        }
        if (!is_constant<T2>::value) {
          adjoint_of(b_arena) += adjB;
        }
      });
}

/**
 * Returns the solution of the system Ax=b when A is triangular and b=I.
 *
 * @tparam T type of elements in A
 * @param A Triangular matrix.
 * @return x = A^-1 .
 * @throws std::domain_error if A is not square
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline var_value<matrix_cl<double>> mdivide_left_tri_low(
    const var_value<T>& A) {
  check_square("mdivide_left_tri_low", "A", A);
  if (A.size() == 0) {
    return A;
  }
  return make_callback_var(
      mdivide_left_tri_low(value_of(A)),
      [A](const vari_value<matrix_cl<double>>& res) {
        matrix_cl<double> res_val_transpose = transpose(res.val());
        matrix_cl<double> adjA
            = res_val_transpose * res.adj() * res_val_transpose;
        adjA.view(matrix_cl_view::Lower);
        adjoint_of(A) -= adjA;
      });
}

}  // namespace math
}  // namespace stan

#endif
#endif
