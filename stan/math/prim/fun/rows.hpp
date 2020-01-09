#ifndef STAN_MATH_PRIM_MAT_FUN_ROWS_HPP
#define STAN_MATH_PRIM_MAT_FUN_ROWS_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return the number of rows in the specified
 * matrix, vector, or row vector.
 *
 * @tparam T type of elements in the matrix
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 *
 * @param[in] m Input matrix, vector, or row vector.
 * @return Number of rows.
 */
template <typename T, int R, int C>
inline int rows(const Eigen::Matrix<T, R, C>& m) {
  return m.rows();
}

}  // namespace math
}  // namespace stan

#endif
