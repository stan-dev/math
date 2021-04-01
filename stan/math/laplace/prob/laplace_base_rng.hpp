#ifndef STAN_MATH_LAPLACE_PROB_LAPLACE_BASE_RNG_HPP
#define STAN_MATH_LAPLACE_PROB_LAPLACE_BASE_RNG_HPP

#include <stan/math/prim/prob/multi_normal_cholesky_rng.hpp>
#include <stan/math/prim/fun/cholesky_decompose.hpp>
#include <stan/math/laplace/laplace_marginal.hpp>

#include <Eigen/Sparse>
#include <Eigen/LU>

namespace stan {
namespace math {

/**
 * In a latent gaussian model,
 *
 *   theta ~ Normal(theta | 0, Sigma(phi, x))
 *   y ~ pi(y | theta, eta)
 *
 * returns a multivariate normal random variate sampled
 * from the gaussian approximation of p(theta_pred | y, phi, x_pred).
 * Note that while the data is observed at x, the new samples
 * are drawn for covariates x_pred.
 * To sample the "original" theta's, set x_pred = x.
 */
template <typename T_theta, typename T_phi, typename T_eta,
          typename T_x, typename T_x_pred,
          typename D, typename K, class RNG>
inline Eigen::VectorXd  // CHECK -- right return type
laplace_base_rng
  (const D& diff_likelihood,
   const K& covariance_function,
   const Eigen::Matrix<T_phi, Eigen::Dynamic, 1>& phi,
   const Eigen::Matrix<T_eta, Eigen::Dynamic, 1>& eta,
   const T_x& x,
   const T_x_pred& x_pred,
   const std::vector<double>& delta,
   const std::vector<int>& delta_int,
   const Eigen::Matrix<T_theta, Eigen::Dynamic, 1>& theta_0,
   RNG& rng,
   std::ostream* msgs = nullptr,
   double tolerance = 1e-6,
   long int max_num_steps = 100,
   int hessian_block_size = 0,
   int compute_W_root = 1) {
  using Eigen::VectorXd;
  using Eigen::MatrixXd;

  VectorXd phi_dbl = value_of(phi);
  VectorXd eta_dbl = value_of(eta);
  Eigen::SparseMatrix<double> W_r;
  MatrixXd L;
  Eigen::PartialPivLU<MatrixXd> LU;
  VectorXd l_grad;
  MatrixXd covariance;
  {
    VectorXd theta;
    VectorXd a;
    double marginal_density
      = laplace_marginal_density(diff_likelihood, covariance_function,
                                 phi_dbl, eta_dbl,
                                 x, delta, delta_int,
                                 covariance, theta, W_r, L, a, l_grad,
                                 LU, value_of(theta_0), msgs,
                                 tolerance, max_num_steps,
                                 hessian_block_size, compute_W_root);
  }

  // Modified R&W method
  MatrixXd covariance_pred = covariance_function(phi_dbl, x_pred,
                                                 delta, delta_int, msgs);
  VectorXd pred_mean = covariance_pred * l_grad;

  Eigen::MatrixXd Sigma;
  if (compute_W_root) {
    Eigen::MatrixXd V_dec = mdivide_left_tri<Eigen::Lower>(L,
                              W_r * covariance_pred);
    Sigma = covariance_pred - V_dec.transpose() * V_dec;
  } else {
    Sigma = covariance_pred
      - covariance_pred * (W_r - W_r * LU.solve(covariance * W_r))
        * covariance_pred;
  }

  return multi_normal_rng(pred_mean, Sigma, rng);
}

}  // namespace math
}  // namespace stan

#endif
