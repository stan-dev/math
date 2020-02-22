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
template <typename EigArr1, typename EigArr2, typename T,
          typename = require_all_eigen_t<EigArr1, EigArr2>>
auto read_cov_matrix(EigArr1&& CPCs, EigArr2&& sds, T& log_prob) {
  return multiply_lower_tri_self_transpose(read_cov_L(CPCs, sds, log_prob));
}

/**
 * Builds a covariance matrix from CPCs and standard deviations
 *
 * @tparam T type of elements in the arrays
 * @param CPCs in (-1, 1)
 * @param sds in (0, inf)
 */
template <typename EigArr1, typename EigArr2,
          typename = require_all_eigen_t<EigArr1, EigArr2>>
auto read_cov_matrix(EigArr1&& CPCs, EigArr2&& sds) {
  using eigen_scalar
      = return_type_t<value_type_t<EigArr1>, value_type_t<EigArr2>>;
  Eigen::Matrix<eigen_scalar, Eigen::Dynamic, Eigen::Dynamic> L
      = sds.matrix().asDiagonal() * read_corr_L(CPCs, sds.rows());
  return multiply_lower_tri_self_transpose(L);
}

}  // namespace math
}  // namespace stan

#endif
