#ifndef STAN_MATH_PRIM_FUN_ROWS_DOT_SELF_HPP
#define STAN_MATH_PRIM_FUN_ROWS_DOT_SELF_HPP

#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Returns the dot product of each row of a matrix with itself.
 *
 * @tparam T type of elements in the matrix
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 *
 * @param x matrix
 */
template <typename T, int R, int C>
inline Eigen::Matrix<T, R, 1> rows_dot_self(const Eigen::Matrix<T, R, C>& x) {
  return x.rowwise().squaredNorm();
}

}  // namespace math
}  // namespace stan

#endif
