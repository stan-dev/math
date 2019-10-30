#ifndef STAN_MATH_OPENCL_DIAGONAL_MULTIPLY_HPP
#define STAN_MATH_OPENCL_DIAGONAL_MULTIPLY_HPP
#ifdef STAN_OPENCL
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/opencl/err/check_opencl.hpp>
#include <stan/math/opencl/kernels/scalar_mul_diagonal.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
namespace stan {
namespace math {
/**
 * Multiplies the diagonal of a matrix on the OpenCL device with the specified
 * scalar.
 *
 * @param A input matrix
 * @param scalar scalar
 * @return copy of the input matrix with the diagonal multiplied by scalar
 */
template <typename T1, typename T2, typename = require_all_arithmetic_t<T1, T2>>
inline matrix_cl<return_type_t<T1, T2>> diagonal_multiply(
    const matrix_cl<T1>& A, const T2 scalar) {
  matrix_cl<return_type_t<T1, T2>> B(A);
  if (B.size() == 0) {
    return B;
  }
  // For rectangular matrices
  int min_dim = B.rows();
  if (B.cols() < min_dim) {
    min_dim = B.cols();
  }
  try {
    opencl_kernels::scalar_mul_diagonal(cl::NDRange(min_dim), B, scalar,
                                        B.rows(), min_dim);
  } catch (const cl::Error& e) {
    check_opencl_error("diagonal_multiply", e);
  }
  return B;
}
}  // namespace math
}  // namespace stan

#endif
#endif
