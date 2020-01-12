#ifndef STAN_MATH_PRIM_FUN_READ_COV_L_HPP
#define STAN_MATH_PRIM_FUN_READ_COV_L_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/constants.hpp>

namespace stan {
namespace math {

/**
 * This is the function that should be called prior to evaluating
 * the density of any elliptical distribution
 *
 * @tparam T type of elements in the arrays
 * @param CPCs on (-1, 1)
 * @param sds on (0, inf)
 * @param log_prob the log probability value to increment with the Jacobian
 * @return Cholesky factor of covariance matrix for specified
 * partial correlations.
 */
template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> read_cov_L(
    const Eigen::Array<T, Eigen::Dynamic, 1>& CPCs,
    const Eigen::Array<T, Eigen::Dynamic, 1>& sds, T& log_prob) {
  size_t K = sds.rows();
  // adjust due to transformation from correlations to covariances
  log_prob += (sum(log(sds)) + LOG_TWO) * K;
  return sds.matrix().asDiagonal() * read_corr_L(CPCs, K, log_prob);
}

}  // namespace math
}  // namespace stan

#endif
