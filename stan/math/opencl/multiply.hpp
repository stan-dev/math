#ifndef STAN_MATH_OPENCL_MULTIPLY_HPP
#define STAN_MATH_OPENCL_MULTIPLY_HPP
#ifdef STAN_OPENCL
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernels/scalar_mul.hpp>
#include <stan/math/opencl/kernels/matrix_multiply.hpp>
#include <Eigen/Dense>

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
          TriangularViewCL triangular_view_B = TriangularViewCL::Entire>
inline auto multiply(const matrix_cl& A, const matrix_cl& B) {
  check_size_match("multiply ((OpenCL))", "A.cols()", A.cols(), "B.rows()",
                   B.rows());
  matrix_cl temp(A.rows(), B.cols());
  if (A.size() == 0 || B.size() == 0) {
    temp.zeros();
    return temp;
  }
  int local = opencl_kernels::matrix_multiply.make_functor.get_opts().at(
      "THREAD_BLOCK_SIZE");
  int Mpad = ((A.rows() + local - 1) / local) * local;
  int Npad = ((B.cols() + local - 1) / local) * local;
  int Kpad = ((A.cols() + local - 1) / local) * local;
  // padding the matrices so the dimensions are divisible with local
  // improves performance and readability because we can omit
  // if statements in the
  // multiply kernel
  matrix_cl tempPad(Mpad, Npad);
  matrix_cl Apad(Mpad, Kpad);
  matrix_cl Bpad(Kpad, Npad);
  opencl_kernels::zeros(cl::NDRange(Mpad, Kpad), Apad.buffer(), Mpad, Kpad,
                        TriangularViewCL::Entire);
  opencl_kernels::zeros(cl::NDRange(Kpad, Npad), Bpad.buffer(), Kpad, Npad,
                        TriangularViewCL::Entire);
  Apad.sub_block<triangular_view_A>(A, 0, 0, 0, 0, A.rows(), A.cols());
  Bpad.sub_block<triangular_view_B>(B, 0, 0, 0, 0, B.rows(), B.cols());
  int wpt = opencl_kernels::matrix_multiply.make_functor.get_opts().at(
      "WORK_PER_THREAD");
  try {
    opencl_kernels::matrix_multiply(
        cl::NDRange(Mpad, Npad / wpt), cl::NDRange(local, local / wpt),
        Apad.buffer(), Bpad.buffer(), tempPad.buffer(), Apad.rows(),
        Bpad.cols(), Bpad.rows(), triangular_view_A, triangular_view_B);
  } catch (cl::Error& e) {
    check_opencl_error("multiply", e);
  }
  // unpadding the result matrix
  temp.sub_block(tempPad, 0, 0, 0, 0, temp.rows(), temp.cols());
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
inline matrix_cl multiply(const matrix_cl& A, const double scalar) {
  matrix_cl temp(A.rows(), A.cols());
  if (A.size() == 0)
    return temp;
  try {
    opencl_kernels::scalar_mul(cl::NDRange(A.rows(), A.cols()), temp.buffer(),
                               A.buffer(), scalar, A.rows(), A.cols());
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
inline auto multiply(const double scalar, const matrix_cl& A) {
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
inline auto multiply(const matrix_cl& A, const matrix_cl& B) {
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
inline matrix_cl operator*(const matrix_cl& A, const matrix_cl& B) {
  return opencl::multiply(A, B);
}
inline matrix_cl operator*(const matrix_cl& B, const double scalar) {
  return multiply(B, scalar);
}
inline matrix_cl operator*(const double scalar, const matrix_cl& B) {
  return multiply(scalar, B);
}
}  // namespace math
}  // namespace stan

#endif
#endif
