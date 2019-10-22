#ifndef STAN_MATH_LAPLACE_LAPLACE_APPROX_POISSON_RNG_HPP
#define STAN_MATH_LAPLACE_LAPLACE_APPROX_POISSON_RNG_HPP

#include <stan/math/laplace/laplace_marginal.hpp>
#include <stan/math/prim/mat/prob/multi_normal_cholesky_rng.hpp>
#include <stan/math/prim/mat/fun/cholesky_decompose.hpp>

namespace stan {
namespace math {

/**
 * In a latent gaussian model,
 * 
 *   theta ~ Normal(0, Sigma(phi))
 *   y ~ log_poisson(theta)
 * 
 * return a multivariate normal random variate sampled
 * from the gaussian approximation of p(theta | y, phi).
 */
template <typename T0, typename T1, class RNG>
inline Eigen::VectorXd  // CHECK -- right return type
laplace_approx_poisson_rng
  (const Eigen::Matrix<T0, Eigen::Dynamic, 1>& theta_0,
   const Eigen::Matrix<T1, Eigen::Dynamic, 1>& phi,
   const std::vector<Eigen::VectorXd>& x,
   const std::vector<int>& n_samples,
   const std::vector<int>& y,
   const Eigen::VectorXd& exposure,
   // const K& covariance_function,
   RNG& rng,
   double tolerance = 1e-6,
   long int max_num_steps = 100) {
  Eigen::MatrixXd covariance;
  Eigen::VectorXd theta;
  Eigen::VectorXd W_root;
  Eigen::MatrixXd L;
  {
    Eigen::VectorXd a;
    Eigen::VectorXd l_grad;
    double marginal_density
      = laplace_marginal_density(theta_0, phi, x,
          diff_poisson_log(to_vector(n_samples), to_vector(y), log(exposure)),
          sqr_exp_kernel_functor(),
          covariance, theta, W_root, L, a, l_grad,
          tolerance, max_num_steps);
  }

  Eigen::MatrixXd V;
  V = mdivide_left_tri<Eigen::Lower>(L,
        diag_pre_multiply(W_root, covariance));

  return multi_normal_cholesky_rng(
    theta,
    cholesky_decompose(covariance - V.transpose() * V),
    rng);
}

}  // namespace math
}  // namespace stan

#endif
