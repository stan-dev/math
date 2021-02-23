#ifndef STAN_MATH_OPENCL_PRIM_REP_MATRIX_HPP
#define STAN_MATH_OPENCL_PRIM_REP_MATRIX_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/err.hpp>

namespace stan {
namespace math {
/** \ingroup opencl
 * Creates a matrix_cl by replicating the given value of arithmetic type.
 *
 * @tparam T type of the result matrix
 * @param x the input value
 * @param n number of rows in the result matrix
 * @param m number of columns in the result matrix
 *
 * @return matrix_cl with replicated value from the input
 *
 * @throw <code>domain_error</code> if the
 * requested dimensions are negative
 *
 */
template <typename T, require_matrix_cl_t<T>* = nullptr>
inline auto rep_matrix(const value_type_t<T>& x, int n, int m) {
  check_nonnegative("rep_matrix (OpenCL)", "rows", n);
  check_nonnegative("rep_matrix (OpenCL)", "cols", m);
  return constant(x, n, m);
}

/** \ingroup opencl
 * Creates a matrix_cl by replicating the input
 * vector or row_vector.  The elements of the
 * vector or row_vector must be of arithmetic type.
 *
 * @tparam T type of elements in the input matrix
 * @param x the input matrix_cl (vector or row_vector)
 * @param m number of rows (if x is a row_vector) or columns
 *  (if x is a vector) in the results matrix
 *
 * @return result matrix with replicated rows or columns
 *
 * @throw <code>domain_error</code> if the
 * requested dimensions are negative
 *
 */
template <typename T, require_arithmetic_t<T>* = nullptr>
inline auto rep_matrix(const matrix_cl<T>& x, int m) {
  check_nonnegative("rep_matrix (OpenCL)", "rows/columns", m);
  check_mat_not_size_one("rep_matrix (OpenCL)", "x", x);
  check_vector("rep_matrix (OpenCL)", "x", x);
  const int N = x.rows() == 1 ? m : x.rows();
  const int M = x.cols() == 1 ? m : x.cols();
  return indexing(x, row_index(N, M) % x.rows(), col_index() % x.cols());
}

}  // namespace math
}  // namespace stan

#endif
#endif
