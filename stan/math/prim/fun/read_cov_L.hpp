#ifndef STAN_MATH_PRIM_FUN_READ_COV_L_HPP
#define STAN_MATH_PRIM_FUN_READ_COV_L_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/read_corr_L.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/constants.hpp>

namespace stan {
namespace math {

/**
 * This is the function that should be called prior to evaluating
 * the density of any elliptical distribution
 *
 * @tparam T_CPCs type of \c T_CPCs vector (must be derived from \c
 * Eigen::ArrayBase and have one compile-time dimension equal to 1)
 * @tparam T_sds type of \c sds vector (must be derived from \c Eigen::ArrayBase
 * and have one compile-time dimension equal to 1)
 * @param CPCs on (-1, 1)
 * @param sds on (0, inf)
 * @param log_prob the log probability value to increment with the Jacobian
 * @return Cholesky factor of covariance matrix for specified
 * partial correlations.
 */
template <typename T_CPCs, typename T_sds,
          require_all_eigen_vector_t<T_CPCs, T_sds>* = nullptr,
          require_vt_same<T_CPCs, T_sds>* = nullptr>
Eigen::Matrix<value_type_t<T_CPCs>, Eigen::Dynamic, Eigen::Dynamic> read_cov_L(
    const T_CPCs& CPCs, const T_sds& sds, value_type_t<T_CPCs>& log_prob) {
  size_t K = sds.rows();
  // adjust due to transformation from correlations to covariances
  log_prob += (sum(log(sds)) + LOG_TWO) * K;
  return make_holder(
      [](const auto& b, const auto& sds) {
        return sds.matrix().asDiagonal() * b;
      },
      read_corr_L(CPCs, K, log_prob), to_ref(sds));
}

}  // namespace math
}  // namespace stan

#endif
