#ifndef STAN_MATH_PRIM_MAT_ERR_IS_ROW_INDEX_HPP
#define STAN_MATH_PRIM_MAT_ERR_IS_ROW_INDEX_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/scal/meta/error_index.hpp>

namespace stan {
namespace math {

/**
 * Return <code>true</code> if the specified index is a valid row of
 * the matrix.
 * This check is 1-indexed by default. This behavior can be changed
 * by setting <code>stan::error_index::value</code>.
 * @tparam T Scalar type, requires class method <code>.rows()</code>
 * @tparam R Compile time rows
 * @tparam C Compile time columns
 * @param y Matrix to test
 * @param i Starting index of row
 * @return <code>true</code> if the index is not out of range
 */
template <typename T_y, int R, int C>
inline bool is_row_index(const Eigen::Matrix<T_y, R, C>& y, size_t i) {
  return i >= error_index::value
         && i < static_cast<size_t>(y.rows()) + error_index::value;
}

}  // namespace math
}  // namespace stan
#endif
