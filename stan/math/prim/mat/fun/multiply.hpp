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

#ifdef STAN_OPENCL
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
template <typename T1, typename T2, enable_if_all_eigen<T1, T2>* = nullptr,
          enable_if_all_eigen_arithmetic<T1, T2>* = nullptr,
          enable_if_not_dot_product<T1, T2>* = nullptr>
inline auto multiply(const T1& m1, const T2& m2) {
  check_multiplicable("multiply", "m1", m1, "m2", m2);
  matrix_cl<double> m1_cl(m1);
  matrix_cl<double> m2_cl(m2);
  matrix_cl<double> m3_cl = m1_cl * m2_cl;
  return from_matrix_cl(m3_cl);
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
template <typename T1, typename T2, enable_if_any_eigen<T1, T2>* = nullptr,
          enable_if_any_not_arithmetic<scalar_type_t<T1>,
                                       scalar_type_t<T2>>* = nullptr,
          enable_if_not_dot_product<T1, T2>* = nullptr,
          enable_if_any_not_contains_var<T1, T2>* = nullptr>
inline auto multiply(const T1& m1, const T2& m2) {
  check_multiplicable("multiply", "m1", m1, "m2", m2);
  return m1 * m2;
}
#else

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
template <typename T1, typename T2,
          enable_if_all_eigen_or_stan_scalar<T1, T2>* = nullptr,
          enable_if_any_not_contains_var<T1, T2>* = nullptr,
          enable_if_not_dot_product<T1, T2>* = nullptr>
inline auto multiply(const T1& m1, const T2& m2) {
  check_multiplicable("multiply", "m1", m1, "m2", m2);
  return m1 * m2;
}
#endif

/**
 * Return the scalar product of the specified row vector and
 * specified column vector.  The return is the same as the dot
 * product. The two vectors must be the same size.
 * @param m1 Row vector.
 * @param m2 Column vector.
 * @return Scalar result of multiplying row vector by column vector.
 * @throw std::domain_error if rv and v are not the same size.
 */
template <typename T1, typename T2,
          enable_if_all_eigen_or_stan_scalar<T1, T2>* = nullptr,
          enable_if_any_not_contains_var<T1, T2>* = nullptr,
          enable_if_dot_product<T1, T2>* = nullptr>
inline auto multiply(const T1& m1, const T2& m2) {
  check_matching_sizes("multiply", "rv", m1, "v", m2);
  return m1.dot(m2);
}

}  // namespace math
}  // namespace stan
#endif
