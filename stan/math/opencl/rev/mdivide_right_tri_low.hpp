#ifndef STAN_MATH_OPENCL_REV_MDIVIDE_RIGHT_TRI_LOW_HPP
#define STAN_MATH_OPENCL_REV_MDIVIDE_RIGHT_TRI_LOW_HPP
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
inline var_value<matrix_cl<double>> mdivide_right_tri_low(T2&& b, T1&& A) {
  check_square("mdivide_right_tri_low", "A", A);
  check_multiplicable("mdivide_right_tri_low", "b", b, "A", A);
  if (A.size() == 0 || b.size() == 0) {
    return var_value<matrix_cl<double>>(matrix_cl<double>(b.rows(), A.cols()));
  }
  arena_t<T1> A_arena = std::forward<T1>(A);
  arena_t<T2> b_arena = std::forward<T2>(b);
  arena_matrix_cl<double> A_tri_inv
      = tri_inverse<matrix_cl_view::Lower>(value_of(A_arena));
  return make_callback_var(
      value_of(b_arena) * A_tri_inv,
      [A_arena, b_arena, A_tri_inv](const vari_value<matrix_cl<double>>& res) {
        matrix_cl<double> adjB = res.adj() * transpose(A_tri_inv);
        if (!is_constant<T1>::value) {
          matrix_cl<double> adjA = transpose(res.val()) * adjB;
          adjA.view(matrix_cl_view::Lower);
          adjoint_of(A_arena) -= adjA;
        }
        if (!is_constant<T2>::value) {
          adjoint_of(b_arena) += adjB;
        }
      });
}

}  // namespace math
}  // namespace stan

#endif
#endif
