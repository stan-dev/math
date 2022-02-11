#ifndef STAN_MATH_REV_FUN_READ_COV_L_HPP
#define STAN_MATH_REV_FUN_READ_COV_L_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/read_cov_L.hpp>
#include <stan/math/rev/fun/read_corr_L.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/constants.hpp>

namespace stan {
namespace math {

/**
 * This is the function that should be called prior to evaluating
 * the density of any elliptical distribution
 *
 * @tparam T_CPCs type of CPCs vector (must be a `var_value<T>` where `T`
 *  inherits from EigenBase)
 * @tparam T_sds type of sds vector (must be a `var_value<T>` where `T`
 *  inherits from EigenBase)
 * @param CPCs on (-1, 1)
 * @param sds on (0, inf)
 * @param log_prob the log probability value to increment with the Jacobian
 * @return Cholesky factor of covariance matrix for specified
 * partial correlations.
 */
template <typename T_CPCs, typename T_sds,
          require_any_var_vector_t<T_CPCs, T_sds>* = nullptr,
          require_vt_same<T_CPCs, T_sds>* = nullptr>
inline auto read_cov_L(const T_CPCs& CPCs, const T_sds& sds,
                       scalar_type_t<T_CPCs>& log_prob) {
  size_t K = sds.rows();
  // adjust due to transformation from correlations to covariances
  log_prob += (sum(log(sds.val())) + LOG_TWO) * K;

  auto corr_L = read_corr_L(CPCs, K, log_prob);
  var_value<Eigen::MatrixXd> res
      = sds.val().matrix().asDiagonal() * corr_L.val();

  reverse_pass_callback([CPCs, sds, corr_L, log_prob, res]() mutable {
    size_t K = sds.size();

    corr_L.adj() += sds.val().matrix().asDiagonal() * res.adj();
    sds.adj() += (res.adj().cwiseProduct(corr_L.val())).rowwise().sum();

    sds.adj() += (K * log_prob.adj() / sds.val().array()).matrix();
  });

  return res;
}

}  // namespace math
}  // namespace stan

#endif
