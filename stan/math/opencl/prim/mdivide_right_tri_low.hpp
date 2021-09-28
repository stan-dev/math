#ifndef STAN_MATH_OPENCL_PRIM_MDIVIDE_RIGHT_TRI_LOW_HPP
#define STAN_MATH_OPENCL_PRIM_MDIVIDE_RIGHT_TRI_LOW_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/err.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/prim/multiply.hpp>
#include <stan/math/opencl/tri_inverse.hpp>

namespace stan {
namespace math {

/**
 * Returns the solution of the system Ax=b where A is a
 * lower triangular matrix.
 *
 * @param A Matrix.
 * @param b Right hand side matrix or vector.
 * @return x = b * tri(A)^-1, solution of the linear system.
 * @throws std::domain_error if A is not square or the rows of b don't
 * match the size of A.
 */
template <typename T1, typename T2,
          require_all_kernel_expressions_t<T1, T2>* = nullptr>
inline matrix_cl<double> mdivide_right_tri_low(const T2& b, const T1& A) {
  check_square("mdivide_right_tri_low (OpenCL)", "A", A);
  check_multiplicable("mdivide_right_tri_low (OpenCL)", "b", b, "A", A);
  if (A.size() == 0 || b.size() == 0) {
    return matrix_cl<double>(b.rows(), A.cols());
  }
  return b * tri_inverse<matrix_cl_view::Lower>(eval(A));
}

}  // namespace math
}  // namespace stan
#endif
#endif
