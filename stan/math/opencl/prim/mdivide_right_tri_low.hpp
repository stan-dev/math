#ifndef STAN_MATH_OPENCL_PRIM_MDIVIDE_RIGHT_TRI_LOW_HPP
#define STAN_MATH_OPENCL_PRIM_MDIVIDE_RIGHT_TRI_LOW_HPP
#ifdef STAN_OPENCL
#include <stan/math/prim/mat/err/check_square.hpp>
#include <stan/math/prim/mat/err/check_multiplicable.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/multiply.hpp>
#include <stan/math/opencl/tri_inverse.hpp>
namespace stan {
namespace math {

/**
 * Returns the solution of the system Ax=b where A is a
 * lower triangular matrix.
 * @param A Matrix.
 * @param b Right hand side matrix or vector.
 * @return x = b * tri(A)^-1, solution of the linear system.
 * @throws std::domain_error if A is not square or the rows of b don't
 * match the size of A.
 */
template <typename T1, typename T2,
          typename = require_all_floating_point_t<T1, T2>>
inline matrix_cl<return_type_t<T1, T2>> mdivide_right_tri_low(
    const matrix_cl<T2>& b, const matrix_cl<T1>& A) {
  check_square("mdivide_right_tri_low (OpenCL)", "A", A);
  check_multiplicable("mdivide_right_tri_low (OpenCL)", "b", b, "A", A);
  return b * tri_inverse<matrix_cl_view::Lower>(A);
}

}  // namespace math
}  // namespace stan
#endif
#endif
