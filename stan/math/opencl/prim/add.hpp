#ifndef STAN_MATH_OPENCL_PRIM_ADD_HPP
#define STAN_MATH_OPENCL_PRIM_ADD_HPP
#ifdef STAN_OPENCL
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/err/check_matching_dims.hpp>
#include <stan/math/opencl/kernels/add.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/prim/meta.hpp>

#include <cl.hpp>

namespace stan {
namespace math {

/**
 * Matrix addition on the OpenCL device
 *
 * @param A first matrix
 * @param B second matrix
 *
 * @return sum of A and B
 *
 * @throw <code>std::invalid_argument</code> if the
 * input matrices do not have matching dimensions
 *
 */
template <typename T1, typename T2, typename = require_all_arithmetic_t<T1, T2>>
inline matrix_cl<return_type_t<T1, T2>> add(const matrix_cl<T1>& A,
                                            const matrix_cl<T2>& B) {
  check_matching_dims("add", "A", A, "B", B);
  matrix_cl<return_type_t<T1, T2>> C(A.rows(), A.cols(),
                                     either(A.view(), B.view()));
  if (C.size() == 0) {
    return C;
  }
  try {
    opencl_kernels::add(cl::NDRange(A.rows(), A.cols()), C, A, B, A.rows(),
                        A.cols(), A.view(), B.view());
  } catch (const cl::Error& e) {
    check_opencl_error("add", e);
  }
  return C;
}

/**
 * Matrix addition on the OpenCL device
 *
 * @param A first matrix
 * @param B second matrix
 *
 * @return sum of A and B
 *
 * @throw <code>std::invalid_argument</code> if the
 * input matrices do not have matching dimensions
 *
 */
template <typename T1, typename T2, typename = require_all_arithmetic_t<T1, T2>>
inline auto operator+(const matrix_cl<T1>& A, const matrix_cl<T2>& B) {
  return add(A, B);
}

}  // namespace math
}  // namespace stan

#endif
#endif
