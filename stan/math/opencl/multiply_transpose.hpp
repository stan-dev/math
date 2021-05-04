#ifndef STAN_MATH_OPENCL_MULTIPLY_TRANSPOSE_HPP
#define STAN_MATH_OPENCL_MULTIPLY_TRANSPOSE_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernels/multiply_transpose.hpp>
#include <stan/math/opencl/err.hpp>

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {
/** \ingroup opencl
 * Computes the product of a square OpenCL matrix with its transpose.
 *
 * Computes the matrix multiplication C = A x A^T
 *
 * @param A input matrix
 * @return the product of the input matrix and its transpose
 *
 */
template <typename T, typename = require_arithmetic_t<T>>
inline matrix_cl<T> multiply_transpose(const matrix_cl<T>& A) {
  matrix_cl<T> temp(A.rows(), A.rows(),
                    A.view() == matrix_cl_view::Diagonal
                        ? matrix_cl_view::Diagonal
                        : matrix_cl_view::Entire);

  if (A.size() == 0) {
    return temp;
  }
  // padding the matrices so the dimensions are divisible with local
  // improves performance because we can omit if statements in the
  // multiply kernel
  int local
      = opencl_kernels::multiply_transpose.get_option("THREAD_BLOCK_SIZE");
  int Mpad = ((A.rows() + local - 1) / local) * local;
  int wpt = opencl_kernels::multiply_transpose.get_option("WORK_PER_THREAD");
  try {
    opencl_kernels::multiply_transpose(cl::NDRange(Mpad, Mpad / wpt),
                                       cl::NDRange(local, local / wpt), A, temp,
                                       A.rows(), A.cols());
  } catch (cl::Error& e) {
    check_opencl_error("multiply self transpose", e);
  }
  return temp;
}
}  // namespace math
}  // namespace stan

#endif
#endif
