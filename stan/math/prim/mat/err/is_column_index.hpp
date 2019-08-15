#ifndef STAN_MATH_PRIM_MAT_ERR_IS_COLUMN_INDEX_HPP
#define STAN_MATH_PRIM_MAT_ERR_IS_COLUMN_INDEX_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

/**
 * Return <code>true</code> no index is invalid column.
 * By default this is a 1-indexed check (as opposed to zero-indexed).
 * Behavior can be changed by setting <code>stan::error_index::value</code>.
 * @tparam T_y Type of scalar, requires class method <code>.cols()</code>
 * @tparam R Number of rows of the matrix
 * @tparam C Number of columns of the matrix
 * @param y Matrix to test
 * @param i Index to check
 * @return <code>true</code> no index is invalid column
 */
template <typename T_y, enable_if_eigen<T_y>* = nullptr>
inline bool is_column_index(const T_y& y, size_t i) {
  return i >= stan::error_index::value
         && i < static_cast<size_t>(y.cols()) + stan::error_index::value;
}

}  // namespace math
}  // namespace stan
#endif
