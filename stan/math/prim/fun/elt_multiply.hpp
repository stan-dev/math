#ifndef STAN_MATH_PRIM_FUN_ELT_MULTIPLY_HPP
#define STAN_MATH_PRIM_FUN_ELT_MULTIPLY_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return the elementwise multiplication of the specified
 * matrices.
 *
 * @tparam Mat1 type of the first matrix or expression
 * @tparam Mat2 type of the second matrix or expression
 *
 * @param m1 First matrix or expression
 * @param m2 Second matrix or expression
 * @return Elementwise product of matrices.
 */
template <typename Mat1, typename Mat2,
          require_all_eigen_t<Mat1, Mat2>* = nullptr,
          require_all_not_st_var<Mat1, Mat2>* = nullptr>
auto elt_multiply(const Mat1& m1, const Mat2& m2) {
  check_matching_dims("elt_multiply", "m1", m1, "m2", m2);
  return m1.cwiseProduct(m2);
}

/**
 * Return the elementwise multiplication of the specified
 * scalars.
 *
 * @tparam Mat1 type of the first scalar
 * @tparam Mat2 type of the second scalar
 *
 * @param a First scalar
 * @param b Second scalar
 * @return Product of scalars
 */
template <typename Scalar1, typename Scalar2,
          require_all_stan_scalar_t<Scalar1, Scalar2>* = nullptr>
auto elt_multiply(const Scalar1& a, const Scalar2& b) {
  return a * b;
}

/**
 * Return specified matrix multiplied by specified scalar.
 *
 * @tparam T1 type of the scalar
 * @tparam T2 type of the matrix or expression
 *
 * @param A scalar
 * @param B matrix
 * @return product of matrix and scalar
 */
template <typename T1, typename T2, require_not_matrix_t<T1>* = nullptr,
          require_matrix_t<T2>* = nullptr,
          require_all_not_st_var<T1, T2>* = nullptr,
          require_not_row_and_col_vector_t<T1, T2>* = nullptr>
inline auto elt_multiply(const T1& A, const T2& B) {
  return (A * B.array()).matrix();
}

/**
 * Return specified matrix multiplied by specified scalar.
 *
 * @tparam T1 type of the matrix or expression
 * @tparam T2 type of the scalar
 *
 * @param A matrix
 * @param B scalar
 * @return product of matrix and scalar
 */
template <typename T1, typename T2, require_matrix_t<T1>* = nullptr,
          require_not_matrix_t<T2>* = nullptr,
          require_all_not_st_var<T1, T2>* = nullptr,
          require_not_row_and_col_vector_t<T1, T2>* = nullptr>
inline auto elt_multiply(const T1& A, const T2& B) {
  return (A.array() * B).matrix();
}



}  // namespace math
}  // namespace stan

#endif
