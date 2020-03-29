#ifndef STAN_MATH_LAPLACE_LAPLACE_APPROX_BERNOULLI_RNG_HPP
#define STAN_MATH_LAPLACE_LAPLACE_APPROX_BERNOULLI_RNG_HPP

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
 * from the gaussian approximation of p(theta | y, phi),
 * where the likelihood is a Bernoulli with logit link.
 */
template <typename K, typename T_phi, typename T_theta,
          typename T_x, class RNG>
inline Eigen::VectorXd  // CHECK -- right return type
  laplace_approx_bernoulli_rng
  (const std::vector<int>& y,
   const std::vector<int>& n_samples,
   const K& covariance_function,
   const Eigen::Matrix<T_phi, Eigen::Dynamic, 1>& phi,
   const T_x x,
   const std::vector<double>& delta,
   const std::vector<int>& delta_int,
   const Eigen::Matrix<T_theta, Eigen::Dynamic, 1>& theta_0,
   RNG& rng,
   std::ostream* msgs = nullptr,
   double tolerance = 1e-6,
   long int max_num_steps = 100) {
    return
    laplace_approx_rng(diff_logistic_log(to_vector(n_samples), to_vector(y)),
                       covariance_function, phi, x, delta, delta_int, theta_0,
                       rng, msgs, tolerance, max_num_steps);
  }

}  // namespace math
}  // namespace stan

#endif
