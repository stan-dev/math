#ifndef STAN_MATH_LAPLACE_PROB_LAPLACE_APPROX_RNG_HPP
#define STAN_MATH_LAPLACE_PROB_LAPLACE_APPROX_RNG_HPP

#include <stan/math/prim/mat/prob/multi_normal_cholesky_rng.hpp>
#include <stan/math/prim/mat/fun/cholesky_decompose.hpp>
#include <stan/math/laplace/laplace_marginal.hpp>

namespace stan {
namespace math {

/**
 * In a latent gaussian model,
 * 
 *   theta ~ Normal(theta | 0, Sigma(phi))
 *   y ~ pi(y | theta)
 * 
 * return a multivariate normal random variate sampled
 * from the gaussian approximation of p(theta | y, phi).
 */
template <typename T0, typename T1, typename D, typename K, class RNG>
inline Eigen::VectorXd  // CHECK -- right return type
laplace_approx_rng
  (const Eigen::Matrix<T0, Eigen::Dynamic, 1>& theta_0,
   const Eigen::Matrix<T1, Eigen::Dynamic, 1>& phi,
   const std::vector<Eigen::VectorXd>& x,
   const D& diff_likelihood,
   const K& covariance_function,
   RNG& rng,
   double tolerance = 1e-6,
   long int max_num_steps = 100) {
  Eigen::VectorXd theta;
  Eigen::VectorXd W_root;
  Eigen::MatrixXd L;
  {
    Eigen::MatrixXd covariance;
    Eigen::VectorXd a;
    Eigen::VectorXd l_grad;
    double marginal_density
      = laplace_marginal_density(value_of(theta_0),
          value_of(phi), x,
          diff_likelihood,
          covariance_function,
          covariance, theta, W_root, L, a, l_grad,
          tolerance, max_num_steps);
  }

  // Modified R&W method
  Eigen::VectorXd W_root_inv = inv(W_root);
  Eigen::MatrixXd V_dec = mdivide_left_tri<Eigen::Lower>(L,
                            diag_matrix(W_root_inv));
  return multi_normal_rng(
    theta,
    diag_matrix(square(W_root_inv)) - - V_dec.transpose() * V_dec,
    rng);
  
  // CHECK -- which method to use? Both seem equivalent.
  // R&W method
  // Eigen::MatrixXd V;
  // V = mdivide_left_tri<Eigen::Lower>(L,
  //       diag_pre_multiply(W_root, covariance));
  // return multi_normal_rng(
  //   theta,
  //   covariance - V.transpose() * V,
  //   rng);
}

}  // namespace math
}  // namespace stan

#endif
