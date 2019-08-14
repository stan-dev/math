#ifndef STAN_MATH_PRIM_MAT_FUN_MULTIPLY_HPP
#define STAN_MATH_PRIM_MAT_FUN_MULTIPLY_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/arr/err/check_matching_sizes.hpp>
#include <stan/math/prim/mat/err/check_multiplicable.hpp>
#ifdef STAN_OPENCL
#include <stan/math/opencl/opencl.hpp>
#endif
#include <type_traits>

namespace stan {
namespace math {

/**
 * Return specified matrix multiplied by specified scalar.
 * @tparam R Row type for matrix.
 * @tparam C Column type for matrix.
 * @param m Matrix.
 * @param c Scalar.
 * @return Product of matrix and scalar.
 */
template <typename T1, typename T2, typename = enable_if_eigen<T1>, typename = enable_if_all_arithmetic<scalar_type_t<T1>, T2>>
inline auto multiply(const T1& m, const T2 c) {
  return c * m;
}

/**
 * Return specified scalar multiplied by specified matrix.
 * @tparam R Row type for matrix.
 * @tparam C Column type for matrix.
 * @param c Scalar.
 * @param m Matrix.
 * @return Product of scalar and matrix.
 */
template <typename T1, typename T2, typename = enable_if_eigen<T2>, typename = enable_if_all_arithmetic<T1, scalar_type_t<T2>>>
inline auto multiply(const T1 c, const T2& m) {
  return c * m;
}

/**
 * Return the product of the specified matrices.  The number of
 * columns in the first matrix must be the same as the number of rows
 * in the second matrix.
 * @param m1 First matrix.
 * @param m2 Second matrix.
 * @return The product of the first and second matrices.
 * @throw std::domain_error if the number of columns of m1 does not match
 *   the number of rows of m2.
 */
template <typename T1, typename T2, typename = enable_if_all_eigen<T1, T2>, typename = enable_if_all_scalar_arithmetic<T1, T2>, typename = enable_if_not_dot_product<T1, T2>>
inline auto multiply(const T1& m1, const T2& m2) {
  check_multiplicable("multiply", "m1", m1, "m2", m2);
#ifdef STAN_OPENCL
  if (m1.rows() * m1.cols() * m2.cols()
      > opencl_context.tuning_opts().multiply_dim_prod_worth_transfer) {
    matrix_cl<double> m1_cl(m1);
    matrix_cl<double> m2_cl(m2);
    matrix_cl<double> m3_cl = m1_cl * m2_cl;
    return from_matrix_cl(m3_cl);
  } else {
    return m1 * m2;
  }
#else
  return m1 * m2;
#endif
}

/**
 * Return the scalar product of the specified row vector and
 * specified column vector.  The return is the same as the dot
 * product.  The two vectors must be the same size.
 * @param rv Row vector.
 * @param v Column vector.
 * @return Scalar result of multiplying row vector by column vector.
 * @throw std::domain_error if rv and v are not the same size.
 */
template <typename T1, typename T2, typename = enable_if_all_scalar_arithmetic<T1, T2>, typename = enable_if_dot_product<T1, T2>> 
inline auto multiply(const T1& rv, const T2& v) {
  check_matching_sizes("multiply", "rv", rv, "v", v);
  return rv.dot(v);
}

/**
 * Return specified matrix multiplied by specified scalar.
 * @tparam R Row type for matrix.
 * @tparam C Column type for matrix.
 * @param m Matrix.
 * @param c Scalar.
 * @return Product of matrix and scalar.
 */
template <typename T1, typename T2, typename = enable_if_all_arithmetic<T1, T2>>
inline return_type_t<T1, T2> multiply(T1 m, T2 c) {
  return c * m;
}

}  // namespace math
}  // namespace stan
#endif
