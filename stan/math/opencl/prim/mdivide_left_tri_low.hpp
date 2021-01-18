#ifndef STAN_MATH_OPENCL_PRIM_MDIVIDE_LEFT_TRI_LOW_HPP
#define STAN_MATH_OPENCL_PRIM_MDIVIDE_LEFT_TRI_LOW_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/err.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/prim/multiply.hpp>
#include <stan/math/opencl/tri_inverse.hpp>

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
template <typename T1, typename T2,
          require_all_kernel_expressions_t<T1, T2>* = nullptr>
inline matrix_cl<double> mdivide_left_tri_low(const T1& A, const T2& b) {
  check_square("mdivide_left_tri_low", "A", A);
  check_multiplicable("mdivide_left_tri_low", "A", A, "b", b);
  if (A.size() == 0 || b.size() == 0) {
    return matrix_cl<double>(A.rows(), b.cols());
  }
  return tri_inverse<matrix_cl_view::Lower>(eval(A)) * b;
}

/**
 * Returns the solution of the system Ax=b when A is triangular and b=I.
 *
 * @tparam T type of elements in A
 * @param A Triangular matrix.
 * @return x = A^-1 .
 * @throws std::domain_error if A is not square
 */
template <typename T, require_all_kernel_expressions_t<T>* = nullptr>
inline matrix_cl<double> mdivide_left_tri_low(const T& A) {
  check_square("mdivide_left_tri_low", "A", A);
  if (A.size() == 0) {
    return A;
  }
  return tri_inverse<matrix_cl_view::Lower>(eval(A));
}

}  // namespace math
}  // namespace stan
#endif
#endif
