#ifndef STAN_MATH_LAPLACE_LAPLACE_APPROX_RNG_HPP
#define STAN_MATH_LAPLACE_LAPLACE_APPROX_RNG_HPP

#include <stan/math/laplace/prob/laplace_approx_rng.hpp>

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
  laplace_approx_poisson_rng
  (const Eigen::Matrix<T0, Eigen::Dynamic, 1>& theta_0,
   const Eigen::Matrix<T1, Eigen::Dynamic, 1>& phi,
   const std::vector<Eigen::VectorXd>& x,
   // const K& covariance_function,
   const std::vector<int>& n_samples,
   const std::vector<int>& y,
   RNG& rng,
   double tolerance = 1e-6,
   long int max_num_steps = 100) {
    return
    laplace_approx_rng(theta_0, phi, x,
                       diff_poisson_log(to_vector(n_samples), to_vector(y)),
                       sqr_exp_kernel_functor(),
                       rng, tolerance, max_num_steps);
}

/**
 * Overload for case where user passes exposure.
 */
template <typename T0, typename T1, typename D, typename K, class RNG>
inline Eigen::VectorXd  // CHECK -- right return type
laplace_approx_poisson_rng
  (const Eigen::Matrix<T0, Eigen::Dynamic, 1>& theta_0,
   const Eigen::Matrix<T1, Eigen::Dynamic, 1>& phi,
   const std::vector<Eigen::VectorXd>& x,
   // const K& covariance_function,
   const std::vector<int>& n_samples,
   const std::vector<int>& y,
   const Eigen::VectorXd& exposure,
   RNG& rng,
   double tolerance = 1e-6,
   long int max_num_steps = 100) {
  return
    laplace_approx_rng(theta_0, phi, x,
                       diff_poisson_log(to_vector(n_samples), to_vector(y),
                                        log(exposure)),
                       sqr_exp_kernel_functor(),
                       rng, tolerance, max_num_steps);
}

}  // namespace math
}  // namespace stan

#endif
