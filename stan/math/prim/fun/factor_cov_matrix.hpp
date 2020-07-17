#ifndef STAN_MATH_PRIM_FUN_FACTOR_COV_MATRIX_HPP
#define STAN_MATH_PRIM_FUN_FACTOR_COV_MATRIX_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/factor_U.hpp>
#include <cstddef>

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
 * @param[out] sds fill this unbounded
 * @return false if any of the diagonals of Sigma are 0
 */
template <typename T_Sigma, typename T_CPCs, typename T_sds,
          require_eigen_t<T_Sigma>* = nullptr,
          require_all_eigen_vector_t<T_CPCs, T_sds>* = nullptr,
          require_all_vt_same<T_Sigma, T_CPCs, T_sds>* = nullptr>
bool factor_cov_matrix(const T_Sigma& Sigma, T_CPCs&& CPCs, T_sds&& sds) {
  using T_scalar = value_type_t<T_Sigma>;
  size_t K = sds.rows();
  const Eigen::Ref<const plain_type_t<T_Sigma>>& Sigma_ref = Sigma;
  sds = Sigma_ref.diagonal().array();
  if ((sds <= 0.0).any()) {
    return false;
  }
  sds = sds.sqrt();

  Eigen::DiagonalMatrix<T_scalar, Eigen::Dynamic> D(K);
  D.diagonal() = sds.inverse();
  sds = sds.log();  // now unbounded

  Eigen::Matrix<T_scalar, Eigen::Dynamic, Eigen::Dynamic> R = D * Sigma_ref * D;
  // to hopefully prevent pivoting due to floating point error
  R.diagonal().setOnes();
  Eigen::LDLT<Eigen::Matrix<T_scalar, Eigen::Dynamic, Eigen::Dynamic>> ldlt;
  ldlt = R.ldlt();
  if (!ldlt.isPositive()) {
    return false;
  }
  Eigen::Matrix<T_scalar, Eigen::Dynamic, Eigen::Dynamic> U = ldlt.matrixU();
  factor_U(U, CPCs);
  return true;
}

}  // namespace math
}  // namespace stan

#endif
