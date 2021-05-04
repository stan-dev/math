#ifndef STAN_MATH_PRIM_PROB_MULTI_GP_LOG_HPP
#define STAN_MATH_PRIM_PROB_MULTI_GP_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/prob/multi_gp_lpdf.hpp>

namespace stan {
namespace math {

/** \ingroup multivar_dists
 * The log of a multivariate Gaussian Process for the given y, Sigma, and
 * w.  y is a dxN matrix, where each column is a different observation and each
 * row is a different output dimension.  The Gaussian Process is assumed to
 * have a scaled kernel matrix with a different scale for each output dimension.
 * This distribution is equivalent to:
 *    for (i in 1:d) row(y, i) ~ multi_normal(0, (1/w[i])*Sigma).
 *
 * @deprecated use <code>multi_gp_lpdf</code>
 *
 * @param y A dxN matrix
 * @param Sigma The NxN kernel matrix
 * @param w A d-dimensional vector of positive inverse scale parameters for each
 * output.
 * @return The log of the multivariate GP density.
 * @throw std::domain_error if Sigma is not square, not symmetric,
 * or not semi-positive definite.
 * @tparam T_y Type of scalar.
 * @tparam T_covar Type of kernel.
 * @tparam T_w Type of weight.
 */
template <bool propto, typename T_y, typename T_covar, typename T_w>
return_type_t<T_y, T_covar, T_w> multi_gp_log(const T_y& y,
                                              const T_covar& Sigma,
                                              const T_w& w) {
  return multi_gp_lpdf<propto>(y, Sigma, w);
}

/** \ingroup multivar_dists
 * @deprecated use <code>multi_gp_lpdf</code>
 */
template <typename T_y, typename T_covar, typename T_w>
inline return_type_t<T_y, T_covar, T_w> multi_gp_log(const T_y& y,
                                                     const T_covar& Sigma,
                                                     const T_w& w) {
  return multi_gp_lpdf<>(y, Sigma, w);
}

}  // namespace math
}  // namespace stan
#endif
