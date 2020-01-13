#ifndef STAN_MATH_PRIM_PROB_MULTI_NORMAL_CHOLESKY_LOG_HPP
#define STAN_MATH_PRIM_PROB_MULTI_NORMAL_CHOLESKY_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/prob/multi_normal_cholesky_lpdf.hpp>

namespace stan {
namespace math {
/** \ingroup multivar_dists
 * The log of the multivariate normal density for the given y, mu, and
 * a Cholesky factor L of the variance matrix.
 * Sigma = LL', a square, semi-positive definite matrix.
 *
 * @deprecated use <code>multi_normal_cholesky_lpdf</code>
 *
 * @param y A scalar vector
 * @param mu The mean vector of the multivariate normal distribution.
 * @param L The Cholesky decomposition of a variance matrix
 * of the multivariate normal distribution
 * @return The log of the multivariate normal density.
 * @throw std::domain_error if LL' is not square, not symmetric,
 * or not semi-positive definite.
 * @tparam T_y Type of scalar.
 * @tparam T_loc Type of location.
 * @tparam T_covar Type of scale.
 */
template <bool propto, typename T_y, typename T_loc, typename T_covar>
return_type_t<T_y, T_loc, T_covar> multi_normal_cholesky_log(const T_y& y,
                                                             const T_loc& mu,
                                                             const T_covar& L) {
  return multi_normal_cholesky_lpdf<propto, T_y, T_loc, T_covar>(y, mu, L);
}

/** \ingroup multivar_dists
 * @deprecated use <code>multi_normal_cholesky_lpdf</code>
 */
template <typename T_y, typename T_loc, typename T_covar>
inline return_type_t<T_y, T_loc, T_covar> multi_normal_cholesky_log(
    const T_y& y, const T_loc& mu, const T_covar& L) {
  return multi_normal_cholesky_lpdf<T_y, T_loc, T_covar>(y, mu, L);
}

}  // namespace math
}  // namespace stan
#endif
