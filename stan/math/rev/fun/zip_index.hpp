#ifndef STAN_MATH_REV_FUN_ZIP_INDEX_HPP
#define STAN_MATH_REV_FUN_ZIP_INDEX_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/to_arena.hpp>
#include <stan/math/prim/err/check_size_match.hpp>
#include <stan/math/prim/fun/zip_index.hpp>
#include <vector>

namespace stan::math {

/**
 * Return a column vector of matrix entries selected by paired one-based row
 * and column indices. Repeated index pairs accumulate adjoints in the same
 * input entry.
 *
 * @tparam EigMat Eigen matrix type inside the input var_value
 * @tparam IdxRows Row index std::vector type with integer elements
 * @tparam IdxCols Column index std::vector type with integer elements
 * @param x Input matrix
 * @param idx_row Row indices
 * @param idx_col Column indices
 * @return A var_value column vector with element i equal to
 * x(idx_row[i] - 1, idx_col[i] - 1), or an empty vector for empty indices
 * @throw std::invalid_argument if the index vectors have different lengths
 * @throw std::out_of_range if an index is outside the matrix bounds
 */
template <typename EigMat, typename IdxRows, typename IdxCols,
          require_eigen_matrix_dynamic_t<EigMat>* = nullptr,
          require_all_std_vector_t<IdxRows, IdxCols>* = nullptr>
inline var_value<Eigen::VectorXd> zip_index(const var_value<EigMat>& x,
                                            IdxRows&& idx_row,
                                            IdxCols&& idx_col) {
  check_size_match("zip_index", "size of idx_row", idx_row.size(),
                   "size of idx_col", idx_col.size());
  if (idx_row.empty()) {
    return var_value<Eigen::VectorXd>(Eigen::VectorXd(0));
  }
  auto rows = to_arena(std::forward<IdxRows>(idx_row));
  auto cols = to_arena(std::forward<IdxCols>(idx_col));
  return make_callback_var(zip_index(x.val(), rows, cols),
                           [x, rows, cols](auto& vi) mutable {
                             zip_index(x.adj(), rows, cols) += vi.adj();
                           });
}

}  // namespace stan::math

#endif
