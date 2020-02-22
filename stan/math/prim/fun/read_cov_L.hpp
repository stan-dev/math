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
 * @tparam Eig1 type for cannonical partial correlations.
 * @tparam Eig2 type for sds
 * @tparam T type for log prob.
 * @param CPCs on (-1, 1)
 * @param sds on (0, inf)
 * @param log_prob the log probability value to increment with the Jacobian
 * @return Cholesky factor of covariance matrix for specified
 * partial correlations.
 */
template <typename Eig1, typename Eig2, typename T,
          typename = require_all_eigen_t<Eig1, Eig2>>
auto read_cov_L(Eig1&& CPCs, Eig2&& sds, T& log_prob) {
  auto K = sds.rows();
  // adjust due to transformation from correlations to covariances
  log_prob += (sum(log(sds)) + LOG_TWO) * K;
  return sds.matrix().asDiagonal() * read_corr_L(CPCs, K, log_prob);
}

}  // namespace math
}  // namespace stan

#endif
