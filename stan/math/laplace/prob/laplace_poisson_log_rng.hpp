#ifndef STAN_MATH_LAPLACE_LAPLACE_APPROX_POISSON_RNG_HPP
#define STAN_MATH_LAPLACE_LAPLACE_APPROX_POISSON_RNG_HPP

#include <stan/math/laplace/prob/laplace_base_rng.hpp>
#include <stan/math/laplace/laplace_likelihood_poisson_log.hpp>
#include <stan/math/laplace/laplace_likelihood_general.hpp>

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
 * where the likelihood is a Poisson with a log link.
 */
template <typename K, typename T0, typename T1, typename T2,  // typename T3,
          typename RNG>
inline Eigen::VectorXd laplace_poisson_log_rng(
    const std::vector<int>& y, const std::vector<int>& n_samples,
    const K& covariance_function,
    const Eigen::Matrix<T0, Eigen::Dynamic, 1>& phi, const T2& x,
    // const T3& x_pred,
    const std::vector<double>& delta, const std::vector<int>& delta_int,
    const Eigen::Matrix<T1, Eigen::Dynamic, 1>& theta_0, RNG& rng,
    std::ostream* msgs = nullptr, double tolerance = 1e-6,
    long int max_num_steps = 100, int hessian_block_size = 0,
    int compute_W_root = 1) {
  Eigen::VectorXd eta_dummy;
  poisson_log_likelihood L;
  return laplace_base_rng(
    diff_likelihood<poisson_log_likelihood>(L, to_vector(y), n_samples, msgs),
    covariance_function, phi, eta_dummy, x, x, delta,
    delta_int, theta_0, rng, msgs, tolerance,
    max_num_steps, hessian_block_size, compute_W_root);
}

/**
 * Overload for case where user passes exposure, ye.
 */
template <typename K, typename T0, typename T1, typename T2,  // typename T3,
          class RNG>
inline Eigen::VectorXd  // CHECK -- right return type
laplace_poisson_log_rng(
    const std::vector<int>& y, const std::vector<int>& n_samples,
    const Eigen::VectorXd& ye, const K& covariance_function,
    const Eigen::Matrix<T0, Eigen::Dynamic, 1>& phi, const T2& x,
    // const T3& x_pred,
    const std::vector<double>& delta, const std::vector<int>& delta_int,
    const Eigen::Matrix<T1, Eigen::Dynamic, 1>& theta_0, RNG& rng,
    std::ostream* msgs = nullptr, double tolerance = 1e-6,
    long int max_num_steps = 100, int hessian_block_size = 0,
    int compute_W_root = 1) {
  Eigen::VectorXd eta_dummy;
  Eigen::VectorXd y_vec = to_vector(y);
  Eigen::VectorXd y_and_ye(y_vec.size() + ye.size());
  y_and_ye << y_vec, ye;
  poisson_log_exposure_likelihood L;
  return laplace_base_rng(
      diff_likelihood<poisson_log_exposure_likelihood>(L, y_and_ye, n_samples,
                                                       msgs),
      covariance_function, phi, eta_dummy, x, x, delta, delta_int, theta_0, rng,
      msgs, tolerance, max_num_steps, hessian_block_size, compute_W_root);
}
}  // namespace math
}  // namespace stan

#endif
