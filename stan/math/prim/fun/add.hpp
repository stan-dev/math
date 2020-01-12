#ifndef STAN_MATH_PRIM_FUN_ADD_HPP
#define STAN_MATH_PRIM_FUN_ADD_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return the sum of the specified matrices.  The two matrices
 * must have the same dimensions.
 *
 * @tparam T1 type of elements in the first matrix
 * @tparam T2 type of elements in the second matrix
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 *
 * @param m1 First matrix.
 * @param m2 Second matrix.
 * @return Sum of the matrices.
 * @throw std::invalid_argument if m1 and m2 do not have the same
 * dimensions.
 */
template <typename T1, typename T2, int R, int C>
inline Eigen::Matrix<return_type_t<T1, T2>, R, C> add(
    const Eigen::Matrix<T1, R, C>& m1, const Eigen::Matrix<T2, R, C>& m2) {
  check_matching_dims("add", "m1", m1, "m2", m2);
  return m1 + m2;
}

/**
 * Return the sum of the specified matrix and specified scalar.
 *
 * @tparam T1 type of elements in the matrix
 * @tparam T2 type of scalar
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 * @param m Matrix.
 * @param c Scalar.
 * @return The matrix plus the scalar.
 */
template <typename T1, typename T2, int R, int C>
inline Eigen::Matrix<return_type_t<T1, T2>, R, C> add(
    const Eigen::Matrix<T1, R, C>& m, const T2& c) {
  return m.array() + c;
}

/**
 * Return the sum of the specified scalar and specified matrix.
 *
 * @tparam T1 type of scalar
 * @tparam T2 type of elements in the matrix
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 * @param c Scalar.
 * @param m Matrix.
 * @return The scalar plus the matrix.
 */
template <typename T1, typename T2, int R, int C>
inline Eigen::Matrix<return_type_t<T1, T2>, R, C> add(
    const T1& c, const Eigen::Matrix<T2, R, C>& m) {
  return c + m.array();
}

}  // namespace math
}  // namespace stan

#endif
