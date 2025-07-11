#ifndef STAN_MATH_PRIM_FUN_SUB_COL_HPP
#define STAN_MATH_PRIM_FUN_SUB_COL_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return a nrows x 1 subcolumn starting at (i-1, j-1).
 *
 * @tparam T type of the matrix
 * @param m Matrix.
 * @param i Starting row + 1.
 * @param j Starting column + 1.
 * @param nrows Number of rows in block.
 * @throw std::out_of_range if either index is out of range.
 */
template <
    typename T, require_matrix_t<T>* = nullptr,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr>
inline auto sub_col(const T& m, size_t i, size_t j, size_t nrows) {
  check_row_index("sub_col", "i", m, i);
  if (nrows > 0) {
    check_row_index("sub_col", "i+nrows-1", m, i + nrows - 1);
  }
  check_column_index("sub_col", "j", m, j);
  return m.col(j - 1).segment(i - 1, nrows);
}

}  // namespace math
}  // namespace stan

#endif
