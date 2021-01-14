#ifndef STAN_MATH_OPENCL_REV_MDIVIDE_LEFT_TRI_LOW_HPP
#define STAN_MATH_OPENCL_REV_MDIVIDE_LEFT_TRI_LOW_HPP
#ifdef STAN_OPENCL

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
  arena_t<T1> A_arena = std::forward<T1>(A);
  arena_t<T2> b_arena = std::forward<T2>(b);
  return make_callback_var(
      mdivide_left_tri_low(value_of(A_arena), value_of(b_arena)),
      [A_arena, b_arena](const vari_value<matrix_cl<double>>& res) {
        matrix_cl<double> adjB = transpose(tri_inverse(A_arena)) * res.adj();
        if (!is_constant<T1>::value) {
          matrix_cl<double> adjA = adjB * transpose(res.adj());
          results(adjoint_of(A_arena), adjoint_of(b_arena))
              += expressions(-adjA, calc_if<!is_constant<T2>::value>(adjB));
        } else {
          adjoint_of(b_arena) += adjB;
        }
      });
}

}  // namespace math
}  // namespace stan

#endif
#endif
