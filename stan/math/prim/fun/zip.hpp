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
template <typename EigMat, typename IdxRows, typename IdxCols,
          require_eigen_matrix_dynamic_t<EigMat>* = nullptr,
          require_all_vector_t<IdxRows, IdxCols>* = nullptr>
inline auto zip(EigMat&& x, IdxRows&& idx_row, IdxCols&& idx_col) {
  check_size_match("zip", "size of idx_row", idx_row.size(), "size of idx_col",
                   idx_col.size());
  return make_holder(
      [](auto&& x_, auto&& idx_row_, auto&& idx_col_) {
        using map_t = Eigen::Map<const Eigen::Array<int, Eigen::Dynamic, 1>>;
        const map_t rows(idx_row_.data(), idx_row_.size());
        const map_t cols(idx_col_.data(), idx_col_.size());
    // If the user turns of range checks do not pay for min and max sweeps
#ifndef STAN_NO_RANGE_CHECKS
        check_range("zip", "minimum row index", x_.rows(), rows.minCoeff());
        check_range("zip", "maximum row index", x_.rows(), rows.maxCoeff());
        check_range("zip", "minimum column index", x_.cols(), cols.minCoeff());
        check_range("zip", "maximum column index", x_.cols(), cols.maxCoeff());
#endif
        const auto linear_idx = (rows.cast<Eigen::Index>() - 1)
                                + (cols.cast<Eigen::Index>() - 1) * x_.rows();
        return x_.reshaped()(linear_idx);
      },
      std::forward<EigMat>(x), std::forward<IdxRows>(idx_row),
      std::forward<IdxCols>(idx_col));
}

}  // namespace math
}  // namespace stan

#endif
