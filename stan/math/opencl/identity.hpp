#ifndef STAN_MATH_OPENCL_IDENTITY_HPP
#define STAN_MATH_OPENCL_IDENTITY_HPP
#ifdef STAN_OPENCL
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/opencl/err/check_opencl.hpp>
#include <stan/math/opencl/kernels/identity.hpp>
#include <stan/math/prim/meta.hpp>
#include <cl.hpp>

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
template <typename T, typename = require_arithmetic_t<T>>
inline matrix_cl<T> identity(int rows_cols) {
  matrix_cl<T> A(rows_cols, rows_cols, matrix_cl_view::Diagonal);
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
