#ifndef STAN_MATH_PRIM_MAT_FUN_COLUMNS_DOT_SELF_HPP
#define STAN_MATH_PRIM_MAT_FUN_COLUMNS_DOT_SELF_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Returns the dot product of each column of a matrix with itself.
 *
 * @tparam T type of elements in the matrix
 * @tparam R number of rows or Eigen::Dynamic
 * @tparam C number of columns or Eigen::Dynamic
 *
 * @param x Matrix.
 */
template <typename T, int R, int C>
inline Eigen::Matrix<T, 1, C> columns_dot_self(
    const Eigen::Matrix<T, R, C>& x) {
  return x.colwise().squaredNorm();
}

}  // namespace math
}  // namespace stan

#endif
