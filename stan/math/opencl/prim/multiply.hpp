#ifndef STAN_MATH_OPENCL_PRIM_MULTIPLY_HPP
#define STAN_MATH_OPENCL_PRIM_MULTIPLY_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/err.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/kernels/matrix_multiply.hpp>
#include <stan/math/opencl/kernels/add.hpp>
#include <stan/math/opencl/scalar_type.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <algorithm>

namespace stan {
namespace math {

/** \ingroup opencl
 * Computes the product of the specified matrices with the option
 * of specifying the triangularity of either input matrices.
 *
 * Computes the matrix multiplication C[M, K] = A[M, N] x B[N, K]
 *
 * @param A first matrix
 * @param B second matrix
 * @tparam partial_view_A specifies whether the matrix A is a
 *  lower/upper triangular or a rectangular matrix
 * @tparam partial_view_B specifies whether the matrix B is a
 *  lower/upper triangular or a rectangular matrix
 * @return the product of the first and second matrix
 *
 * @throw <code>std::invalid_argument</code> if the
 *   number of columns in A and rows in B do not match
 */

template <typename T1, typename T2,
          typename = require_all_kernel_expressions_and_none_scalar_t<T1, T2>>
inline matrix_cl<return_type_t<T1, T2>> multiply(const T1& A, const T2& B) {
  check_size_match("multiply ((OpenCL))", "A.cols()", A.cols(), "B.rows()",
                   B.rows());
  if (A.size() == 0 || B.size() == 0) {
    return constant(0.0, A.rows(), B.cols());
  }
  matrix_cl<return_type_t<T1, T2>> temp(A.rows(), B.cols(),
                                        either(A.view(), B.view()));
  if (A.rows() == 1) {
    const int local_size
        = opencl_kernels::row_vector_matrix_multiply.get_option("LOCAL_SIZE_");
    try {
      opencl_kernels::row_vector_matrix_multiply(
          cl::NDRange(temp.cols() * local_size), cl::NDRange(local_size),
          A.eval(), B.eval(), temp, B.rows(), B.cols(), A.view(), B.view());
    } catch (cl::Error& e) {
      check_opencl_error("row_vector - matrix multiply", e);
    }
    return temp;
  }
  if (B.cols() == 1) {
    temp = matrix_vector_multiply(A, B);
    return temp;
  }
  int local = opencl_kernels::matrix_multiply.get_option("THREAD_BLOCK_SIZE");
  const int Mpad = ((A.rows() + local - 1) / local) * local;
  const int Npad = ((B.cols() + local - 1) / local) * local;
  const int wpt = opencl_kernels::matrix_multiply.get_option("WORK_PER_THREAD");
  const int wgs = Mpad / local * Npad / local;
  const int split = std::min(
      A.cols() / local,
      (opencl_context.tuning_opts().multiply_wgs_per_compute_unit
           * static_cast<int>(opencl_context.device()[0]
                                  .getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>())
       + wgs - 1)
          / wgs);
  try {
    if (split <= 1) {
      opencl_kernels::matrix_multiply(cl::NDRange(Mpad, Npad / wpt),
                                      cl::NDRange(local, local / wpt), A.eval(),
                                      B.eval(), temp, A.rows(), B.cols(),
                                      B.rows(), A.view(), B.view());
    } else {
      matrix_cl<return_type_t<T1, T2>> tempSplit(A.rows(), B.cols() * split);
      opencl_kernels::matrix_multiply(cl::NDRange(Mpad, Npad / wpt, split),
                                      cl::NDRange(local, local / wpt, 1),
                                      A.eval(), B.eval(), tempSplit, A.rows(),
                                      B.cols(), B.rows(), A.view(), B.view());
      opencl_kernels::add_batch(cl::NDRange(A.rows(), B.cols()), temp,
                                tempSplit, A.rows(), B.cols(), split);
    }
  } catch (cl::Error& e) {
    check_opencl_error("multiply", e);
  }
  return temp;
}

/**
 * Matrix multiplication of two kernel generator expressions. Evaluates both
 * expressions before calculating the matrix product.
 * @tparam T_a type of first expression
 * @tparam T_b type of second expression
 * @param a first expression
 * @param b second expression
 * @return Matrix product of given arguments
 */
template <typename T_a, typename T_b,
          typename = require_all_kernel_expressions_and_none_scalar_t<T_a, T_b>>
inline matrix_cl<return_type_t<T_a, T_b>> operator*(const T_a& a,
                                                    const T_b& b) {
  // no need for perfect forwarding as operations are evaluated
  return multiply(as_operation_cl(a).eval(), as_operation_cl(b).eval());
}

/**
 * Matrix multiplication of a scalar and a  kernel generator expressions.
 * @tparam T_a type of scalar
 * @tparam T_b type of the kernel generator expression
 * @param a scalar
 * @param b expression
 * @return Matrix product of given arguments
 */
template <typename T_a, typename T_b, require_stan_scalar_t<T_a>* = nullptr,
          require_all_kernel_expressions_and_none_scalar_t<T_b>* = nullptr,
          require_all_not_var_t<T_a, T_b>* = nullptr>
inline matrix_cl<return_type_t<T_a, T_b>> multiply(const T_a& a, const T_b& b) {
  return a * b;
}

/**
 * Matrix multiplication of a scalar and a  kernel generator expressions.
 * @tparam T_a type of the kernel generator expression
 * @tparam T_b type of scalar
 * @param a expression
 * @param b scalar
 * @return Matrix product of given arguments
 */
template <typename T_a, typename T_b, require_stan_scalar_t<T_b>* = nullptr,
          require_all_kernel_expressions_and_none_scalar_t<T_a>* = nullptr,
          require_all_not_var_t<T_a, T_b>* = nullptr>
inline matrix_cl<return_type_t<T_a, T_b>> multiply(const T_a& a, const T_b& b) {
  return a * b;
}

}  // namespace math
}  // namespace stan
#endif
#endif
