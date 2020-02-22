#ifndef STAN_MATH_PRIM_FUN_FACTOR_COV_MATRIX_HPP
#define STAN_MATH_PRIM_FUN_FACTOR_COV_MATRIX_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/factor_U.hpp>
#include <cstddef>
#include <tuple>

namespace stan {
namespace math {

/**
 * This function is intended to make starting values, given a
 * covariance matrix Sigma
 *
 * The transformations are hard coded as log for standard
 * deviations and Fisher transformations (atanh()) of CPCs
 *
 * @tparam T type of elements in the matrix and arrays
 * @param[in] Sigma covariance matrix
 * @param[out] CPCs fill this unbounded (does not resize)
 * @param[out] sds fill this unbounded (does not resize)
 * @return false if any of the diagonals of Sigma are 0
 */
template <typename EigMat, typename = require_eigen_t<EigMat>>
inline auto factor_cov_matrix(EigMat&& Sigma) {
  using eigen_scalar = value_type_t<EigMat>;
  Eigen::Array<eigen_scalar, Eigen::Dynamic, 1> sds = Sigma.diagonal().array();
  size_t K = sds.rows();
  if ((sds <= 0.0).any()) {
    throw_domain_error("factor_cov_matrix", "failed on y", "", "");
  }
  sds = sds.sqrt();
  Eigen::DiagonalMatrix<eigen_scalar, Eigen::Dynamic> D(K);
  D.diagonal() = sds.inverse();
  sds = sds.log();  // now unbounded

  Eigen::Matrix<eigen_scalar, Eigen::Dynamic, Eigen::Dynamic> R = D * Sigma * D;
  // to hopefully prevent pivoting due to floating point error
  R.diagonal().setOnes();
  auto ldlt = R.ldlt();
  if (!ldlt.isPositive()) {
    throw_domain_error("factor_cov_matrix", "failed on y", "", "");
  }
  Eigen::Array<eigen_scalar, Eigen::Dynamic, 1> CPCs = factor_U(ldlt.matrixU());
  return std::make_tuple(CPCs, sds);
}

}  // namespace math
}  // namespace stan

#endif
