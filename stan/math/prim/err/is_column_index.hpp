#ifndef STAN_MATH_PRIM_ERR_IS_COLUMN_INDEX_HPP
#define STAN_MATH_PRIM_ERR_IS_COLUMN_INDEX_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return <code>true</code> if column index is in bounds.
 * By default this is a 1-indexed check (as opposed to zero-indexed).
 * Behavior can be changed by setting <code>stan::error_index::value</code>.
 * @tparam T_y Type of scalar, requires class method <code>.cols()</code>
 * @tparam R number of rows or Eigen::Dynamic
 * @tparam C number of columns or Eigen::Dynamic
 * @param y matrix to test
 * @param i column index to check
 * @return <code>true</code> if column index is in bounds
 */
template <typename T_y, int R, int C>
inline bool is_column_index(const Eigen::Matrix<T_y, R, C>& y, size_t i) {
  return i >= stan::error_index::value
         && i < static_cast<size_t>(y.cols()) + stan::error_index::value;
}

}  // namespace math
}  // namespace stan
#endif
