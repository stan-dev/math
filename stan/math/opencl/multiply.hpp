#ifndef STAN_MATH_OPENCL_MULTIPLY_HPP
#define STAN_MATH_OPENCL_MULTIPLY_HPP
#ifdef STAN_OPENCL
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/err/check_opencl.hpp>
#include <stan/math/opencl/kernels/scalar_mul.hpp>
#include <stan/math/opencl/kernels/matrix_multiply.hpp>
#include <stan/math/opencl/kernels/add.hpp>
#include <stan/math/opencl/sub_block.hpp>
#include <stan/math/opencl/zeros.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {
namespace opencl {
/**
 * Computes the product of the specified matrices with the option
 * of specifying the triangularity of either input matrices.
 *
 * Computes the matrix multiplication C[M, K] = A[M, N] x B[N, K]
 *
 * @param A first matrix
 * @param B second matrix
 * @tparam triangular_view_A specifies whether the matrix A is a
 *  lower/upper triangular or a rectangular matrix
 * @tparam triangular_view_B specifies whether the matrix B is a
 *  lower/upper triangular or a rectangular matrix
 * @return the product of the first and second matrix
 *
 * @throw <code>std::invalid_argument</code> if the
 *   number of columns in A and rows in B do not match
 */

template <TriangularViewCL triangular_view_A = TriangularViewCL::Entire,
          TriangularViewCL triangular_view_B = TriangularViewCL::Entire,
          typename T1, typename T2, typename = enable_if_all_arithmetic<T1, T2>>
inline matrix_cl<return_type_t<T1, T2>> multiply(const matrix_cl<T1>& A,
                                                 const matrix_cl<T2>& B) {
  check_size_match("multiply ((OpenCL))", "A.cols()", A.cols(), "B.rows()",
                   B.rows());
  matrix_cl<return_type_t<T1, T2>> temp(A.rows(), B.cols());
  if (A.size() == 0 || B.size() == 0) {
    temp.zeros();
    return temp;
  }
  if (A.rows() == 1) {
    const int local_size
        = opencl_kernels::row_vector_matrix_multiply.make_functor.get_opts().at(
            "LOCAL_SIZE_");
    try {
      opencl_kernels::row_vector_matrix_multiply(
          cl::NDRange(temp.cols() * local_size), cl::NDRange(local_size), A, B,
          temp, B.rows(), B.cols(), triangular_view_A, triangular_view_B);
    } catch (cl::Error& e) {
      check_opencl_error("row_vector - matrix multiply", e);
    }
    return temp;
  }
  if (B.cols() == 1) {
    try {
      opencl_kernels::matrix_vector_multiply(
          cl::NDRange(temp.rows()), A, B, temp, A.rows(), A.cols(),
          triangular_view_A, triangular_view_B);
    } catch (cl::Error& e) {
      check_opencl_error("matrix - vector multiply", e);
    }
    return temp;
  }
  int local = opencl_kernels::matrix_multiply.make_functor.get_opts().at(
      "THREAD_BLOCK_SIZE");
  const int Mpad = ((A.rows() + local - 1) / local) * local;
  const int Npad = ((B.cols() + local - 1) / local) * local;
  const int wpt = opencl_kernels::matrix_multiply.make_functor.get_opts().at(
      "WORK_PER_THREAD");
  int split = A.cols() / std::sqrt(A.rows() * B.cols());
  if (split > 20) {
    split = 20;
  }
  // when there result matrix is large, there is no benefit of splitting
  // as the number of created threads is large enough to occupy all
  // compute units in the OpenCL device
  if (temp.size() > opencl_context.tuning_opts().multiply_split_upper_limit) {
    split = 1;
  }
  try {
    if (split <= 1) {
      opencl_kernels::matrix_multiply(cl::NDRange(Mpad, Npad / wpt),
                                      cl::NDRange(local, local / wpt), A, B,
                                      temp, A.rows(), B.cols(), B.rows(),
                                      triangular_view_A, triangular_view_B);
    } else {
      matrix_cl<return_type_t<T1, T2>> tempSplit(A.rows(), B.cols() * split);
      opencl_kernels::matrix_multiply(cl::NDRange(Mpad, Npad / wpt, split),
                                      cl::NDRange(local, local / wpt, 1), A, B,
                                      tempSplit, A.rows(), B.cols(), B.rows(),
                                      triangular_view_A, triangular_view_B);
      opencl_kernels::add_batch(cl::NDRange(A.rows(), B.cols()), temp,
                                tempSplit, A.rows(), B.cols(), split);
    }
  } catch (cl::Error& e) {
    check_opencl_error("multiply", e);
  }
  return temp;
}
}  // namespace opencl

/**
 * Multiplies the specified matrix on the OpenCL device
 * with the specified scalar.
 *
 * @param A matrix
 * @param scalar scalar
 * @return matrix multipled with scalar
 */
template <typename T1, typename T2, typename = enable_if_all_arithmetic<T1, T2>>
inline matrix_cl<return_type_t<T1, T2>> multiply(const matrix_cl<T1>& A,
                                                 const T2 scalar) {
  matrix_cl<return_type_t<T1, T2>> temp(A.rows(), A.cols());
  if (A.size() == 0)
    return temp;
  try {
    opencl_kernels::scalar_mul(cl::NDRange(A.rows(), A.cols()), temp, A, scalar,
                               A.rows(), A.cols());
  } catch (const cl::Error& e) {
    check_opencl_error("multiply scalar", e);
  }
  return temp;
}

/**
 * Multiplies the specified matrix on the OpenCL device
 * with the specified scalar.
 *
 * @param scalar scalar
 * @param A matrix
 * @return matrix multipled with scalar
 */
template <typename T1, typename T2, typename = enable_if_all_arithmetic<T1, T2>>
inline matrix_cl<return_type_t<T1, T2>> multiply(const T1 scalar,
                                                 const matrix_cl<T2>& A) {
  return multiply(A, scalar);
}

/**
 * Computes the product of the specified matrices.
 *
 * Computes the matrix multiplication C[M, K] = A[M, N] x B[N, K]
 *
 * @param A first matrix
 * @param B second matrix
 * @return the product of the first and second matrix
 *
 * @throw <code>std::invalid_argument</code> if the
 *   number of columns in A and rows in B do not match
 */
template <typename T1, typename T2, typename = enable_if_all_arithmetic<T1, T2>>
inline matrix_cl<return_type_t<T1, T2>> multiply(const matrix_cl<T1>& A,
                                                 const matrix_cl<T2>& B) {
  return opencl::multiply(A, B);
}

/**
 * Templated product operator for OpenCL matrices.
 *
 * Computes the matrix multiplication C[M, K] = A[M, N] x B[N, K].
 *
 * @param A A matrix or scalar
 * @param B A matrix or scalar
 * @return the product of the first and second arguments
 *
 * @throw <code>std::invalid_argument</code> if the
 *   number of columns in A and rows in B do not match
 */
template <typename T1, typename T2, typename = enable_if_all_arithmetic<T1, T2>>
inline matrix_cl<return_type_t<T1, T2>> operator*(const matrix_cl<T1>& A,
                                                  const matrix_cl<T2>& B) {
  return opencl::multiply(A, B);
}
template <typename T1, typename T2, typename = enable_if_all_arithmetic<T1, T2>>
inline matrix_cl<return_type_t<T1, T2>> operator*(const matrix_cl<T1>& B,
                                                  const T2 scalar) {
  return multiply(B, scalar);
}
template <typename T1, typename T2, typename = enable_if_all_arithmetic<T1, T2>>
inline matrix_cl<return_type_t<T1, T2>> operator*(const T1 scalar,
                                                  const matrix_cl<T2>& B) {
  return multiply(scalar, B);
}
}  // namespace math
}  // namespace stan

#endif
#endif
