#ifndef STAN_MATH_PRIM_FUN_BLOCK_HPP
#define STAN_MATH_PRIM_FUN_BLOCK_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return a nrows x ncols submatrix starting at (i-1, j-1).
 *
 * @tparam T type of elements in the matrix
 * @param m Matrix.
 * @param i Starting row.
 * @param j Starting column.
 * @param nrows Number of rows in block.
 * @param ncols Number of columns in block.
 * @throw std::out_of_range if either index is out of range.
 */
template <typename T, require_matrix_t<T>* = nullptr>
inline auto block(const T& m, size_t i, size_t j, size_t nrows, size_t ncols) {
  check_row_index("block", "i", m, i);
  check_row_index("block", "i+nrows-1", m, i + nrows - 1);
  check_column_index("block", "j", m, j);
  check_column_index("block", "j+ncols-1", m, j + ncols - 1);
  return m.block(i - 1, j - 1, nrows, ncols);
}

}  // namespace math
}  // namespace stan

#endif
