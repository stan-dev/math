#ifndef STAN_MATH_OPENCL_DIVIDE_COLUMNS_HPP
#define STAN_MATH_OPENCL_DIVIDE_COLUMNS_HPP
#ifdef STAN_OPENCL
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/prim/scal/err/check_size_match.hpp>
#include <stan/math/opencl/kernels/divide_columns.hpp>
#include <CL/cl.hpp>

namespace stan {
namespace math {

/**
 * Divides matrix columns element-wise by a vector
 *
 * @param A first matrix
 * @param B Vector of elements to divide each column of @c A by element-wise.
 *
 * @return element-wise division of @c A by @c B.
 *
 * @throw <code>std::invalid_argument</code> if the
 * input matrices do not have matching dimensions
 *
 */
inline void divide_columns(const matrix_cl& A, const matrix_cl& B) try {
  check_size_match("divide_columns", "A", A.rows(), "B", B.size());
  check_size_match("divide_columns", "A", B.cols(), "B", 1);
  if (A.size() == 0 || B.size() == 0) {
    return;
  }
  opencl_kernels::divide_columns_vec(cl::NDRange(A.size()), A, B, B.size());
  return;
  } catch (const cl::Error& e) {
    check_opencl_error("divide_columns", e);
    return;
  }


/**
* Takes matrix @c A element-wise
 *
 * @param A Matrix to divide
 * @param divisor scalar to divide by
 *
 * @return element-wise division of @c A by @c divisor.
 *
 */
inline void divide_columns(const matrix_cl& A, const double& divisor) try {
  if (A.size() == 0) {
    return;
  }
  opencl_kernels::divide_columns_scalar(cl::NDRange(A.size()), A, divisor);
  return;
  } catch (const cl::Error& e) {
    check_opencl_error("divide_columns", e);
    return;
  }


}  // namespace math
}  // namespace stan

#endif
#endif
