#ifndef STAN_MATH_OPENCL_PRIM_SUBTRACT_HPP
#define STAN_MATH_OPENCL_PRIM_SUBTRACT_HPP
#ifdef STAN_OPENCL
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/opencl/kernels/subtract.hpp>
#include <stan/math/opencl/err/check_matching_dims.hpp>
#include <stan/math/opencl/err/check_opencl.hpp>
#include <stan/math/prim/meta.hpp>
#include <cl.hpp>

namespace stan {
namespace math {

/**
 * Matrix subtraction on the OpenCL device
 * Subtracts the second matrix
 * from the first matrix and stores
 * the result in the third matrix (C=A-B)
 *
 * @param A first matrix
 * @param B second matrix
 *
 * @return subtraction result matrix
 *
 * @throw <code>std::invalid_argument</code> if the
 * input matrices do not have matching dimensions.
 *
 */
template <typename T1, typename T2, typename = require_all_arithmetic_t<T1, T2>>
inline matrix_cl<return_type_t<T1, T2>> subtract(const matrix_cl<T1>& A,
                                                 const matrix_cl<T2>& B) {
  check_matching_dims("subtract ((OpenCL))", "A", A, "B", B);
  matrix_cl<return_type_t<T1, T2>> C(A.rows(), A.cols(),
                                     either(A.view(), B.view()));
  if (A.size() == 0) {
    return C;
  }
  try {
    opencl_kernels::subtract(cl::NDRange(A.rows(), A.cols()), C, A, B, A.rows(),
                             A.cols(), A.view(), B.view());
  } catch (cl::Error& e) {
    check_opencl_error("subtract", e);
  }
  return C;
}

/**
 * Matrix subtraction on the OpenCL device
 * Subtracts the second matrix
 * from the first matrix and stores
 * the result in the third matrix (C=A-B)
 *
 * @param A first matrix
 * @param B second matrix
 *
 * @return subtraction result matrix
 *
 * @throw <code>std::invalid_argument</code> if the
 * input matrices do not have matching dimensions.
 *
 */
template <typename T1, typename T2, typename = require_all_arithmetic_t<T1, T2>>
inline auto operator-(const matrix_cl<T1>& A, const matrix_cl<T2>& B) {
  return subtract(A, B);
}
}  // namespace math
}  // namespace stan

#endif
#endif
