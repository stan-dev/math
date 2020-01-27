#ifndef STAN_MATH_PRIM_FUN_MULTIPLY_HPP
#define STAN_MATH_PRIM_FUN_MULTIPLY_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/dot_product.hpp>
#ifdef STAN_OPENCL
#include <stan/math/opencl/opencl.hpp>
#endif
#include <type_traits>

namespace stan {
namespace math {

/**
 * Return specified matrix multiplied by specified scalar.
 *
 * @tparam T1 type of the matrix or expression
 * @tparam T2 type of the scalar
 *
 * @param m matrix
 * @param c scalar
 * @return product of matrix and scalar
 */
template <typename T1, typename T2, typename = require_eigen_t<T1>,
          typename = require_stan_scalar_t<T2>>
inline auto multiply(const T1& m, T2 c) {
  return (c * m).eval();
}

/**
 * Return specified scalar multiplied by specified matrix.
 *
 * @tparam T1 type of the scalar
 * @tparam T2 type of the matrix or expression
 *
 * @param c scalar
 * @param m matrix
 * @return product of scalar and matrix
 */
template <typename T1, typename T2, typename = require_stan_scalar_t<T1>,
          typename = require_eigen_t<T2>>
inline auto multiply(T1 c, const T2& m) {
  return (c * m).eval();
}

/**
 * Return the product of the specified matrices. The number of
 * columns in the first matrix must be the same as the number of rows
 * in the second matrix.
 *
 * @tparam T1 type of the first matrix or expression
 * @tparam T2 type of the second matrix or expression
 *
 * @param m1 first matrix or expression
 * @param m2 second matrix or expression
 * @return the product of the first and second matrices
 * @throw <code>std::invalid_argument</code> if the number of columns of m1 does
 * not match the number of rows of m2.
 */
template <typename T1, typename T2,
          typename = require_all_eigen_vt<std::is_arithmetic, T1, T2>,
          typename
          = require_any_not_same_t<double, value_type_t<T1>, value_type_t<T2>>,
          typename = require_not_eigen_row_and_col_t<T1, T2>>
inline auto multiply(const T1& m1, const T2& m2) {
  check_multiplicable("multiply", "m1", m1, "m2", m2);
  return (m1 * m2).eval();
}

/**
 * Return the product of the specified matrices. The number of
 * columns in the first matrix must be the same as the number of rows
 * in the second matrix. If scalar of matrices is \c double OpenCL
 * implementation can be used.
 *
 * @tparam T1 type of the first matrix or expression
 * @tparam T2 type of the second matrix or expression
 *
 * @param m1 first matrix or expression
 * @param m2 second matrix or expression
 * @return the product of the first and second matrices
 * @throw <code>std::invalid_argument</code> if the number of columns of m1 does
 * not match the number of rows of m2.
 */
template <typename T1, typename T2, typename = require_all_eigen_t<T1, T2>,
          typename
          = require_all_same_t<double, value_type_t<T1>, value_type_t<T2>>,
          typename = require_not_eigen_row_and_col_t<T1, T2>>
inline auto multiply(const T1& m1, const T2& m2) -> decltype((m1 * m2).eval()) {
  check_multiplicable("multiply", "m1", m1, "m2", m2);
#ifdef STAN_OPENCL
  if (m1.rows() * m1.cols() * m2.cols()
      > opencl_context.tuning_opts().multiply_dim_prod_worth_transfer) {
    matrix_cl<double> m1_cl(m1);
    matrix_cl<double> m2_cl(m2);
    matrix_cl<double> m3_cl = m1_cl * m2_cl;
    return from_matrix_cl<T1::RowsAtCompileTime, T2::ColsAtCompileTime>(m3_cl);
  } else {
    return (m1 * m2).eval();
  }
#else
  return (m1 * m2).eval();
#endif
}

/**
 * Return the scalar product of the specified row vector and
 * specified column vector.  The return is the same as the dot
 * product.  The two vectors must be the same size.
 *
 * @tparam T1 type of the row vector
 * @tparam T2 type of the column vector
 *
 * @param rv row vector
 * @param v column vector
 * @return scalar result of multiplying row vector by column vector
 * @throw <code>std::invalid_argument</code> if rv and v are not the same size
 */
template <typename T1, typename T2,
          typename
          = require_all_not_var_t<scalar_type_t<T1>, scalar_type_t<T2>>,
          typename = require_eigen_row_and_col_t<T1, T2>>
inline auto multiply(const T1& rv, const T2& v) {
  check_matching_sizes("multiply", "rv", rv, "v", v);
  return dot_product(rv, v);
}

/**
 * Return product of scalars.
 *
 * @tparam T1 type of first scalar
 * @tparam T2 type of second scalar
 * @param m scalar
 * @param c scalar
 * @return product
 */
template <typename T1, typename T2,
          typename = require_all_stan_scalar_t<T1, T2>>
inline return_type_t<T1, T2> multiply(T1 m, T2 c) {
  return c * m;
}

}  // namespace math
}  // namespace stan

#endif
