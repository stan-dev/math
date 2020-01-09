#ifndef STAN_MATH_PRIM_MAT_FUN_SUBTRACT_HPP
#define STAN_MATH_PRIM_MAT_FUN_SUBTRACT_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return the result of subtracting the second specified matrix
 * from the first specified matrix.  The return scalar type is the
 * promotion of the input types.
 *
 * @tparam T1 type of elements in the first matrix
 * @tparam T2 type of elements in the second matrix
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 *
 * @param m1 First matrix.
 * @param m2 Second matrix.
 * @return Difference between first matrix and second matrix.
 */
template <typename T1, typename T2, int R, int C>
inline Eigen::Matrix<return_type_t<T1, T2>, R, C> subtract(
    const Eigen::Matrix<T1, R, C>& m1, const Eigen::Matrix<T2, R, C>& m2) {
  check_matching_dims("subtract", "m1", m1, "m2", m2);
  return m1 - m2;
}

template <typename T1, typename T2, int R, int C>
inline Eigen::Matrix<return_type_t<T1, T2>, R, C> subtract(
    const T1& c, const Eigen::Matrix<T2, R, C>& m) {
  return c - m.array();
}

template <typename T1, typename T2, int R, int C>
inline Eigen::Matrix<return_type_t<T1, T2>, R, C> subtract(
    const Eigen::Matrix<T1, R, C>& m, const T2& c) {
  return m.array() - c;
}

}  // namespace math
}  // namespace stan

#endif
