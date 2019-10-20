#ifndef STAN_MATH_OPENCL_PRIM_MDIVIDE_LEFT_TRI_LOW_HPP
#define STAN_MATH_OPENCL_PRIM_MDIVIDE_LEFT_TRI_LOW_HPP
#ifdef STAN_OPENCL
#include <stan/math/prim/mat/err/check_square.hpp>
#include <stan/math/prim/mat/err/check_multiplicable.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/multiply.hpp>
#include <stan/math/opencl/tri_inverse.hpp>
namespace stan {
namespace math {

/**
 * Returns the solution of the system Ax=b when A is lower triangular.
 * @tparam T1 type of elements in A
 * @tparam T2 type of elements in b
 * @param A Triangular matrix.
 * @param b Right hand side matrix or vector.
 * @return x = A^-1 b, solution of the linear system.
 * @throws std::domain_error if A is not square or the rows of b don't
 * match the size of A.
 */
template <typename T1, typename T2,
          typename = require_all_floating_point_t<T1, T2>>
inline matrix_cl<return_type_t<T1, T2>> mdivide_left_tri_low(
    const matrix_cl<T1>& A, const matrix_cl<T2>& b) {
  check_square("mdivide_left_tri_low", "A", A);
  check_multiplicable("mdivide_left_tri_low", "A", A, "b", b);
  return tri_inverse<matrix_cl_view::Lower>(A) * b;
}

/**
 * Returns the solution of the system Ax=b when A is triangular and b=I.
 * @tparam T type of elements in A
 * @tparam R1 number of rows in A
 * @tparam C1 number of columns in A
 * @param A Triangular matrix.
 * @return x = A^-1 .
 * @throws std::domain_error if A is not square
 */
template <typename T, typename = require_all_floating_point_t<T>>
inline matrix_cl<T> mdivide_left_tri_low(const matrix_cl<T>& A) {
  check_square("mdivide_left_tri_low", "A", A);
  return tri_inverse<matrix_cl_view::Lower>(A);
}

}  // namespace math
}  // namespace stan
#endif
#endif
