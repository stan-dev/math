#ifndef STAN_MATH_OPENCL_IDENTITY_HPP
#define STAN_MATH_OPENCL_IDENTITY_HPP
#ifdef STAN_OPENCL
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/err/check_opencl.hpp>
#include <stan/math/opencl/kernels/identity.hpp>
#include <CL/cl.hpp>

namespace stan {
namespace math {

/**
 * Returns the identity matrix stored on the OpenCL device
 *
 * @param rows_cols the number of rows and columns
 *
 * @return the identity matrix
 *
 */
template <typename T, typename std::enable_if_t<std::is_arithmetic<T>::value, int> = 0>
inline matrix_cl<T> identity(int rows_cols) {
  matrix_cl<T> A(rows_cols, rows_cols);
  if (rows_cols == 0) {
    return A;
  }
  try {
    opencl_kernels::identity(cl::NDRange(A.rows(), A.cols()), A, A.rows(),
                             A.cols());
  } catch (const cl::Error& e) {
    check_opencl_error("identity", e);
  }
  return A;
}
}  // namespace math
}  // namespace stan

#endif
#endif
