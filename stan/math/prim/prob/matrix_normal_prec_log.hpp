#ifndef STAN_MATH_PRIM_PROB_MATRIX_NORMAL_PREC_LOG_HPP
#define STAN_MATH_PRIM_PROB_MATRIX_NORMAL_PREC_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/prob/matrix_normal_prec_lpdf.hpp>

namespace stan {
namespace math {
/** \ingroup multivar_dists
 * The log of the matrix normal density for the given y, mu, Sigma and D
 * where Sigma and D are given as precision matrices, not covariance
 * matrices.
 *
 * @deprecated use <code>matrix_normal_prec_lpdf</code>
 *
 * @param y An mxn matrix.
 * @param Mu The mean matrix.
 * @param Sigma The mxm inverse covariance matrix (i.e., the precision
 *   matrix) of the rows of y.
 * @param D The nxn inverse covariance matrix (i.e., the precision
 *   matrix) of the columns of y.
 * @return The log of the matrix normal density.
 * @throw std::domain_error if Sigma or D are not square, not symmetric,
 * or not semi-positive definite.
 * @tparam T_y Type of scalar.
 * @tparam T_Mu Type of location.
 * @tparam T_Sigma Type of Sigma.
 * @tparam T_D Type of D.
 */
template <bool propto, typename T_y, typename T_Mu, typename T_Sigma,
          typename T_D>
return_type_t<T_y, T_Mu, T_Sigma, T_D> matrix_normal_prec_log(
    const Eigen::Matrix<T_y, Eigen::Dynamic, Eigen::Dynamic>& y,
    const Eigen::Matrix<T_Mu, Eigen::Dynamic, Eigen::Dynamic>& Mu,
    const Eigen::Matrix<T_Sigma, Eigen::Dynamic, Eigen::Dynamic>& Sigma,
    const Eigen::Matrix<T_D, Eigen::Dynamic, Eigen::Dynamic>& D) {
  return matrix_normal_prec_lpdf<propto, T_y, T_Mu, T_Sigma, T_D>(y, Mu, Sigma,
                                                                  D);
}

/** \ingroup multivar_dists
 * @deprecated use <code>matrix_normal_prec_lpdf</code>
 */
template <typename T_y, typename T_Mu, typename T_Sigma, typename T_D>
return_type_t<T_y, T_Mu, T_Sigma, T_D> matrix_normal_prec_log(
    const Eigen::Matrix<T_y, Eigen::Dynamic, Eigen::Dynamic>& y,
    const Eigen::Matrix<T_Mu, Eigen::Dynamic, Eigen::Dynamic>& Mu,
    const Eigen::Matrix<T_Sigma, Eigen::Dynamic, Eigen::Dynamic>& Sigma,
    const Eigen::Matrix<T_D, Eigen::Dynamic, Eigen::Dynamic>& D) {
  return matrix_normal_prec_lpdf<T_y, T_Mu, T_Sigma, T_D>(y, Mu, Sigma, D);
}

}  // namespace math
}  // namespace stan
#endif
