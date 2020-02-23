#ifndef STAN_MATH_PRIM_FUN_COV_MATRIX_FREE_LKJ_HPP
#define STAN_MATH_PRIM_FUN_COV_MATRIX_FREE_LKJ_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return the vector of unconstrained partial correlations and
 * deviations that transform to the specified covariance matrix.
 *
 * <p>The constraining transform is defined as for
 * <code>cov_matrix_constrain(Matrix, size_t)</code>.  The
 * inverse first factors out the deviations, then applies the
 * freeing transform of <code>corr_matrix_free(Matrix&)</code>.
 *
 * @tparam T type of elements in the matrix
 * @param y Covariance matrix to free.
 * @return Vector of unconstrained values that transforms to the
 * specified covariance matrix.
 * @throw std::domain_error if the correlation matrix has no
 *    elements or is not a square matrix.
 * @throw std::runtime_error if the correlation matrix cannot be
 *    factorized by factor_cov_matrix()
 */
template <typename EigMat, typename = require_eigen_t<EigMat>>
inline auto cov_matrix_free_lkj(EigMat&& y) {
  using Eigen::Dynamic;
  using eigen_scalar = value_type_t<EigMat>;
  using size_type = index_type_t<Eigen::Matrix<eigen_scalar, Dynamic, Dynamic>>;

  check_nonzero_size("cov_matrix_free_lkj", "y", y);
  check_square("cov_matrix_free_lkj", "y", y);
  size_type k = y.rows();
  size_type k_choose_2 = (k * (k - 1)) / 2;
  auto ret_vals = factor_cov_matrix(std::forward<EigMat>(y));
  Eigen::Matrix<eigen_scalar, Dynamic, 1> x(k_choose_2 + k);
  x.head(k_choose_2) = std::get<0>(ret_vals).head(k_choose_2);
  x.segment(k_choose_2, k) = std::get<1>(ret_vals).head(k);
  return x;
}

}  // namespace math
}  // namespace stan

#endif
