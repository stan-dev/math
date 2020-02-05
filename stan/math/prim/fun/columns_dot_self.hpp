#ifndef STAN_MATH_PRIM_FUN_COLUMNS_DOT_SELF_HPP
#define STAN_MATH_PRIM_FUN_COLUMNS_DOT_SELF_HPP

#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Returns the dot product of each column of a matrix with itself.
 *
 * @tparam T type of elements in the matrix
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 *
 * @param x Matrix.
 * @return Vector containing the dot product of each column of the matrix
 * with itself.
 */
template <typename T, int R, int C>
inline Eigen::Matrix<T, 1, C> columns_dot_self(
    const Eigen::Matrix<T, R, C>& x) {
  return x.colwise().squaredNorm();
}

}  // namespace math
}  // namespace stan

#endif
