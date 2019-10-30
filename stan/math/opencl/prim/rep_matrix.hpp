#ifndef STAN_MATH_OPENCL_PRIM_REP_MATRIX_HPP
#define STAN_MATH_OPENCL_PRIM_REP_MATRIX_HPP
#ifdef STAN_OPENCL
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/opencl/kernels/rep_matrix.hpp>
#include <stan/math/opencl/err/check_opencl.hpp>
#include <stan/math/prim/scal/err/check_nonnegative.hpp>
#include <stan/math/prim/scal/err/check_greater.hpp>
#include <stan/math/prim/scal/err/domain_error.hpp>
#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {
/**
 * Creates a matrix_cl by replicating the value of
 * the only element in the input 1x1 matrix_cl. The
 * element must be of arithmetic type.
 *
 * @tparam T type of elements in the input matrix
 * @param x the input 1x1 matrix_cl
 * @param n number of rows in the results matrix
 * @param m number of columns in the results matrix
 *
 * @return matrix_cl with replicated value from the input matrix
 *
 * @throw <code>domain_error</code> if the
 * requested dimensions are negative
 * @throw <code>invalid_argument</code> if input
 * element is not a matrix_cl of size 1
 *
 */
template <typename T, typename = require_arithmetic_t<T>>
inline matrix_cl<T> rep_matrix(const matrix_cl<T>& x, int n, int m) {
  check_nonnegative("rep_matrix (OpenCL)", "rows", n);
  check_nonnegative("rep_matrix (OpenCL)", "cols", m);
  check_mat_size_one("rep_matrix (OpenCL)", "x", x);
  matrix_cl<T> A(n, m);
  if (A.size() == 0) {
    return A;
  }
  try {
    opencl_kernels::rep_matrix(cl::NDRange(A.rows(), A.cols()), A, x, A.rows(),
                               A.cols(), x.rows(), x.cols(), A.view());
  } catch (cl::Error& e) {
    check_opencl_error("rep_matrix", e);
  }
  return A;
}

/**
 * Creates a matrix_cl by replicating the input
 * vector or row_vector.  The elements of the
 * vector or row_vector must be of arithmetic type.
 *
 * @tparam T type of elements in the input matrix
 * @param x the input matrix_cl (vector or row_vector)
 * @param m number of rows (if x is a row_vecotr) or columns
 *  (if x is a vector) in the results matrix
 *
 * @return result matrix with replicated rows or columns
 *
 * @throw <code>domain_error</code> if the
 * requested dimensions are negative
 *
 */
template <typename T, typename = require_arithmetic_t<T>>
inline matrix_cl<T> rep_matrix(const matrix_cl<T>& x, int m) {
  check_nonnegative("rep_matrix (OpenCL)", "rows/columns", m);
  check_mat_not_size_one("rep_matrix (OpenCL)", "x", x);
  check_vector("rep_matrix (OpenCL)", "x", x);
  const int N = x.rows() == 1 ? m : x.rows();
  const int M = x.cols() == 1 ? m : x.cols();
  matrix_cl<T> A(N, M);
  if (A.size() == 0) {
    return A;
  }
  try {
    opencl_kernels::rep_matrix(cl::NDRange(A.rows(), A.cols()), A, x, A.rows(),
                               A.cols(), x.rows(), x.cols(), A.view());
  } catch (cl::Error& e) {
    check_opencl_error("rep_matrix", e);
  }
  return A;
}

}  // namespace math
}  // namespace stan

#endif
#endif
