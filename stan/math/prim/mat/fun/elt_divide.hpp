#ifndef STAN_MATH_PRIM_MAT_FUN_ELT_DIVIDE_HPP
#define STAN_MATH_PRIM_MAT_FUN_ELT_DIVIDE_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return the elementwise division of the specified matrices.
 *
 * @tparam T1 type of elements in first matrix
 * @tparam T2 type of elements in second matrix
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 *
 * @param m1 First matrix
 * @param m2 Second matrix
 * @return Elementwise division of matrices.
 */
template <typename T1, typename T2, int R, int C>
Eigen::Matrix<return_type_t<T1, T2>, R, C> elt_divide(
    const Eigen::Matrix<T1, R, C>& m1, const Eigen::Matrix<T2, R, C>& m2) {
  check_matching_dims("elt_divide", "m1", m1, "m2", m2);

  return m1.array() / m2.array();
}

/**
 * Return the elementwise division of the specified matrix
 * by the specified scalar.
 *
 * @tparam T1 type of elements in the matrix
 * @tparam T2 type of the scalar
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 *
 * @param m matrix
 * @param s scalar
 * @return Elementwise division of a scalar by matrix.
 */
template <typename T1, typename T2, int R, int C>
Eigen::Matrix<return_type_t<T1, T2>, R, C> elt_divide(
    const Eigen::Matrix<T1, R, C>& m, T2 s) {
  return m.array() / s;
}

/**
 * Return the elementwise division of the specified scalar
 * by the specified matrix.
 *
 * @tparam T1 type of the scalar
 * @tparam T2 type of elements in the matrix
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 *
 * @param s scalar
 * @param m matrix
 * @return Elementwise division of a scalar by matrix.
 */
template <typename T1, typename T2, int R, int C>
Eigen::Matrix<return_type_t<T1, T2>, R, C> elt_divide(
    T1 s, const Eigen::Matrix<T2, R, C>& m) {
  return s / m.array();
}

}  // namespace math
}  // namespace stan

#endif
