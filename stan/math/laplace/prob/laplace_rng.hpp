#ifndef STAN_MATH_LAPLACE_LAPLACE_RNG_HPP
#define STAN_MATH_LAPLACE_LAPLACE_RNG_HPP

#include <stan/math/laplace/prob/laplace_rng.hpp>

namespace stan {
namespace math {

/**
 * In a latent gaussian model,
 *
 *   theta ~ Normal(theta | 0, Sigma(phi))
 *   y ~ pi(y | theta)
 *
 * return a multivariate normal random variate sampled
 * from the gaussian approximation of p(theta | y, phi)
 * where the log likelihood is given by L_f.
 */
template <typename T_phi, typename T_eta, typename T_theta, typename T_x,
          typename K, typename L, typename RNG>
inline Eigen::VectorXd
  laplace_rng
  (const L& L_f,
   const Eigen::Matrix<T_eta, Eigen::Dynamic, 1>& eta,
   const Eigen::VectorXd& delta_L,
   const std::vector<int>& delta_int_L,
   const K& K_f,
   const Eigen::Matrix<T_phi, Eigen::Dynamic, 1>& phi,
   const T_x& x,
   const std::vector<double>& delta_K,
   const std::vector<int>& delta_int_K,
   const Eigen::Matrix<T_theta, Eigen::Dynamic, 1>& theta_0,
   RNG& rng,
   std::ostream* msgs_L = nullptr,
   std::ostream* msgs_K = nullptr,
   double tolerance = 1e-6,
   long int max_num_steps = 100,
   int hessian_block_size = 0,
   int compute_W_root = 1) {
     return
       laplace_base_rng(
         diff_likelihood<L>(L_f, delta_L, delta_int_L, msgs_L),
         K_f, phi, eta,
         x, x, delta_K, delta_int_K, theta_0,
         rng, msgs_K, tolerance, max_num_steps,
         hessian_block_size, compute_W_root);
  }

}  // namespace math
}  // namespace stan

#endif
