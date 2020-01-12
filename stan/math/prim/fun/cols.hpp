#ifndef STAN_MATH_PRIM_FUN_COLS_HPP
#define STAN_MATH_PRIM_FUN_COLS_HPP

#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return the number of columns in the specified
 * matrix, vector, or row vector.
 *
 * @tparam T type of elements in the matrix
 * @tparam R number of rows or Eigen::Dynamic
 * @tparam C number of columns or Eigen::Dynamic
 * @param[in] m Input matrix, vector, or row vector.
 * @return Number of columns.
 */
template <typename T, int R, int C>
inline int cols(const Eigen::Matrix<T, R, C>& m) {
  return m.cols();
}

}  // namespace math
}  // namespace stan

#endif
