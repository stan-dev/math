#ifndef STAN_MATH_PRIM_FUN_READ_COV_MATRIX_HPP
#define STAN_MATH_PRIM_FUN_READ_COV_MATRIX_HPP

#include <stan/math/prim/fun/read_cov_L.hpp>
#include <stan/math/prim/fun/multiply_lower_tri_self_transpose.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * A generally worse alternative to call prior to evaluating the
 * density of an elliptical distribution
 *
 * @tparam T type of elements in the arrays
 * @param CPCs on (-1, 1)
 * @param sds on (0, inf)
 * @param log_prob the log probability value to increment with the Jacobian
 * @return Covariance matrix for specified partial correlations.
 */
template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> read_cov_matrix(
    const Eigen::Array<T, Eigen::Dynamic, 1>& CPCs,
    const Eigen::Array<T, Eigen::Dynamic, 1>& sds, T& log_prob) {
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> L
      = read_cov_L(CPCs, sds, log_prob);
  return multiply_lower_tri_self_transpose(L);
}

/**
 *
 * Builds a covariance matrix from CPCs and standard deviations
 *
 * @tparam T type of elements in the arrays
 * @param CPCs in (-1, 1)
 * @param sds in (0, inf)
 */
template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> read_cov_matrix(
    const Eigen::Array<T, Eigen::Dynamic, 1>& CPCs,
    const Eigen::Array<T, Eigen::Dynamic, 1>& sds) {
  size_t K = sds.rows();
  Eigen::DiagonalMatrix<T, Eigen::Dynamic> D(K);
  D.diagonal() = sds;
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> L = D * read_corr_L(CPCs, K);
  return multiply_lower_tri_self_transpose(L);
}

}  // namespace math
}  // namespace stan

#endif
