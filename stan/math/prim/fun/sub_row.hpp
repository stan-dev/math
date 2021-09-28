#ifndef STAN_MATH_PRIM_FUN_SUB_ROW_HPP
#define STAN_MATH_PRIM_FUN_SUB_ROW_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return a 1 x ncols subrow starting at (i-1, j-1).
 *
 * @tparam T type of the matrix
 * @param m Matrix Input matrix.
 * @param i Starting row + 1.
 * @param j Starting column + 1.
 * @param ncols Number of columns in block.
 * @throw std::out_of_range if either index is out of range.
 */
template <
    typename T, require_matrix_t<T>* = nullptr,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr>
inline auto sub_row(const T& m, size_t i, size_t j, size_t ncols) {
  check_row_index("sub_row", "i", m, i);
  check_column_index("sub_row", "j", m, j);
  if (ncols > 0) {
    check_column_index("sub_col", "j+ncols-1", m, j + ncols - 1);
  }
  return m.row(i - 1).segment(j - 1, ncols);
}

}  // namespace math
}  // namespace stan

#endif
