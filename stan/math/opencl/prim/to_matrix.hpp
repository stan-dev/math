#ifndef STAN_MATH_OPENCL_PRIM_TO_MATRIX_HPP
#define STAN_MATH_OPENCL_PRIM_TO_MATRIX_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernel_generator.hpp>

namespace stan {
namespace math {

/** \ingroup opencl
 * Returns input matrix.
 *
 * @tparam T_x type of the matrix
 *
 * @param x matrix
 * @return the matrix representation of the input
 */
template <typename T_x,
          require_nonscalar_prim_or_rev_kernel_expression_t<T_x>* = nullptr>
inline T_x to_matrix(T_x&& x) {
  return std::forward<T_x>(x);
}

/**
 * Returns a matrix representation of a vector or matrix in column-major
 * order with the specified number of rows and columns.
 *
 * @tparam T_x type of the matrix
 *
 * @param x matrix
 * @param m rows
 * @param n columns
 * @return Reshaped input matrix
 * @throw <code>std::invalid_argument</code> if the sizes
 * do not match
 */
template <typename T_x,
          require_all_kernel_expressions_and_none_scalar_t<T_x>* = nullptr>
inline matrix_cl<return_type_t<T_x>> to_matrix(const T_x& x, int m, int n) {
  using res_scal = return_type_t<T_x>;
  check_size_match("to_matrix", "rows * columns", "", m * n, "input size", "",
                   x.size());
  matrix_cl<res_scal> res(m, n);
  matrix_cl<res_scal> tmp(res.buffer(), x.rows(), x.cols());
  tmp = x;
  for (cl::Event e : tmp.write_events()) {
    res.add_write_event(e);
  }
  return res;
}

/**
 * Returns a matrix representation of the vector or matrix in column-major or
 * row major order with the specified number of rows and columns.
 *
 * @tparam T_x type of the matrix
 * @param x matrix
 * @param m rows
 * @param n columns
 * @param col_major column-major indicator
 * if 1, output matrix is transversed in column-major order
 * if 0, output matrix is transversed in row-major order
 * @return the matrix representation of the input
 * @throw <code>std::invalid_argument</code>
 * if the sizes do not match
 */
template <typename T_x,
          require_nonscalar_prim_or_rev_kernel_expression_t<T_x>* = nullptr>
inline auto to_matrix(const T_x& x, int m, int n, bool col_major)
    -> decltype(to_matrix(x, m, n)) {
  if (col_major) {
    return to_matrix(x, m, n);
  } else {
    return transpose(to_matrix(transpose(x), n, m));
  }
}

}  // namespace math
}  // namespace stan
#endif
#endif
