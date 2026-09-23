#ifndef STAN_MATH_PRIM_FUN_ZIP_HPP
#define STAN_MATH_PRIM_FUN_ZIP_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Return a vector of matrix values obtained by zipping
 * two N-dimensional integer vectors to form tuples of
 * corresponding elements.
 *
 * @tparam EigMat type of the matrix
 * @param x input matrix
 * @param idx_row row index vector
 * @param idx_col column index vector
 * @return column vector with element i equal to x(idx_row[i]-1, idx_col[i]-1)
 * @throw std::invalid_argument if the index vectors are of different lengths
 * @throw std::out_of_range if any index is out of the valid range of the matrix
 */
template <typename EigMat, require_eigen_matrix_dynamic_t<EigMat>* = nullptr>
inline auto zip(const EigMat& x, const std::vector<int>& idx_row,
                const std::vector<int>& idx_col) {
  check_size_match("zip", "size of idx_row", idx_row.size(), "size of idx_col",
                   idx_col.size());

  using map_t = Eigen::Map<const Eigen::Array<int, Eigen::Dynamic, 1>>;

  const map_t rows(idx_row.data(), idx_row.size());
  const map_t cols(idx_col.data(), idx_col.size());

  check_range("zip", "minimum row index", x.rows(), rows.minCoeff());
  check_range("zip", "maximum row index", x.rows(), rows.maxCoeff());
  check_range("zip", "minimum column index", x.cols(), cols.minCoeff());
  check_range("zip", "maximum column index", x.cols(), cols.maxCoeff());

  const auto linear_idx = (rows.cast<Eigen::Index>() - 1)
                          + (cols.cast<Eigen::Index>() - 1) * x.rows();

  return x.reshaped()(linear_idx);
}

}  // namespace math
}  // namespace stan

#endif
