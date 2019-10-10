#ifndef STAN_MATH_OPENCL_PRIM_MULTIPLY_HPP
#define STAN_MATH_OPENCL_PRIM_MULTIPLY_HPP
#ifdef STAN_OPENCL
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/err/check_opencl.hpp>
#include <stan/math/opencl/kernels/scalar_mul.hpp>
#include <stan/math/opencl/kernels/matrix_multiply.hpp>
#include <stan/math/opencl/kernels/add.hpp>
#include <stan/math/opencl/multiply.hpp>
#include <stan/math/opencl/sub_block.hpp>
#include <stan/math/opencl/zeros.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {
/**
 * Multiplies the specified matrix on the OpenCL device
 * with the specified scalar.
 *
 * @tparam T1 type of elements in the matrix
 * @tparam T2 type of scalar
 * @param A matrix
 * @param scalar scalar
 * @return matrix multipled with scalar
 */
template <typename T1, typename T2, typename = require_all_arithmetic_t<T1, T2>>
inline matrix_cl<return_type_t<T1, T2>> multiply(const matrix_cl<T1>& A,
                                                 const T2 scalar) {
  matrix_cl<return_type_t<T1, T2>> temp(A.rows(), A.cols(), A.view());
  if (A.size() == 0) {
    return temp;
  }
  try {
    opencl_kernels::scalar_mul(cl::NDRange(A.rows(), A.cols()), temp, A, scalar,
                               A.rows(), A.cols(), A.view());
  } catch (const cl::Error& e) {
    check_opencl_error("multiply scalar", e);
  }
  return temp;
}

/**
 * Multiplies the specified matrix on the OpenCL device
 * with the specified scalar.
 *
 * @tparam T1 type of scalar
 * @tparam T2 type of elements in the matrix
 * @param scalar scalar
 * @param A matrix
 * @return matrix multipled with scalar
 */
template <typename T1, typename T2, typename = require_all_arithmetic_t<T1, T2>>
inline matrix_cl<return_type_t<T1, T2>> multiply(const T1 scalar,
                                                 const matrix_cl<T2>& A) {
  return multiply(A, scalar);
}

/**
 * Computes the product of the specified matrices.
 *
 * Computes the matrix multiplication C[M, K] = A[M, N] x B[N, K]
 *
 * @tparam T1 type of elements in matrix A
 * @tparam T2 type of elements in matrix B
 * @param A first matrix
 * @param B second matrix
 * @return the product of the first and second matrix
 *
 * @throw <code>std::invalid_argument</code> if the
 *   number of columns in A and rows in B do not match
 */
template <typename T1, typename T2, typename = require_all_arithmetic_t<T1, T2>>
inline matrix_cl<return_type_t<T1, T2>> multiply(const matrix_cl<T1>& A,
                                                 const matrix_cl<T2>& B) {
  return opencl::multiply(A, B);
}

/**
 * Templated product operator for a OpenCL matrices.
 *
 * Computes the matrix multiplication C[M, K] = A[M, N] x B[N, K].
 * @tparam T1 type of elements in matrix A
 * @tparam T2 type of elements in matrix B
 * @param A A matrix or scalar
 * @param B A matrix or scalar
 * @return the product of the first and second arguments
 *
 * @throw <code>std::invalid_argument</code> if the
 *   number of columns in A and rows in B do not match
 */
template <typename T1, typename T2, typename = require_all_arithmetic_t<T1, T2>>
inline matrix_cl<return_type_t<T1, T2>> operator*(const matrix_cl<T1>& A,
                                                  const matrix_cl<T2>& B) {
  return opencl::multiply(A, B);
}
/**
 * Templated product operator for an OpenCL matrix and a scalar.
 *
 * @tparam T1 type of elements in the matrix
 * @tparam T2 type of scalar
 * @param B input matrix
 * @param scalar input scalar
 * @return the product of the first and second arguments
 *
 * @throw <code>std::invalid_argument</code> if the
 *   number of columns in A and rows in B do not match
 */
template <typename T1, typename T2, typename = require_all_arithmetic_t<T1, T2>>
inline matrix_cl<return_type_t<T1, T2>> operator*(const matrix_cl<T1>& B,
                                                  const T2 scalar) {
  return multiply(B, scalar);
}
/**
 * Templated product operator for an OpenCL matrix and a scalar.
 *
 * @tparam T1 type of scalar
 * @tparam T2 type of elements in the matrix
 * @param scalar input scalar
 * @param B input matrix
 * @return the product of the first and second arguments
 *
 * @throw <code>std::invalid_argument</code> if the
 *   number of columns in A and rows in B do not match
 */
template <typename T1, typename T2, typename = require_all_arithmetic_t<T1, T2>>
inline matrix_cl<return_type_t<T1, T2>> operator*(const T1 scalar,
                                                  const matrix_cl<T2>& B) {
  return multiply(scalar, B);
}
}  // namespace math
}  // namespace stan

#endif
#endif
