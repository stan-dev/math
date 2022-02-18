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
template <typename CovarFun, typename T1, class RNG, typename TupleData,
          typename... Args>
inline Eigen::VectorXd laplace_poisson_log_rng(
    const std::vector<int>& y, const std::vector<int>& n_samples,
    CovarFun&& covariance_function,
    const Eigen::Matrix<T1, Eigen::Dynamic, 1>& theta_0,
    const TupleData& data_tuple, RNG& rng, std::ostream* msgs = nullptr,
    const double tolerance = 1e-6, const long int max_num_steps = 100,
    const int hessian_block_size = 0, const int solver = 1,
    const int max_steps_line_search = 0,
    Args&&... args) {
  Eigen::VectorXd eta_dummy;
  poisson_log_likelihood L;
  return laplace_base_rng(
      diff_likelihood<poisson_log_likelihood>(L, to_vector(y), n_samples, msgs),
      covariance_function, eta_dummy, theta_0, data_tuple, rng, msgs, tolerance,
      max_num_steps, hessian_block_size, solver,
      max_steps_line_search, std::forward<Args>(args)...);
}

/**
 * Overload for case where user passes exposure, ye.
 */
template <typename CovarFun, typename T1, class RNG, typename... TrainData,
          typename TupleData, typename... Args>
inline Eigen::VectorXd  // CHECK -- right return type
laplace_poisson_log_rng(
    const std::vector<int>& y, const std::vector<int>& n_samples,
    const Eigen::VectorXd& ye, CovarFun&& covariance_function,
    const Eigen::Matrix<T1, Eigen::Dynamic, 1>& theta_0,
    const TupleData& data_tuple, RNG& rng, std::ostream* msgs = nullptr,
    const double tolerance = 1e-6, const long int max_num_steps = 100,
    const int hessian_block_size = 0, const int solver = 1,
    const int max_steps_line_search = 0,
    Args&&... args) {
  Eigen::VectorXd eta_dummy;
  Eigen::VectorXd y_vec = to_vector(y);
  Eigen::VectorXd y_and_ye(y_vec.size() + ye.size());
  y_and_ye << y_vec, ye;
  poisson_log_exposure_likelihood L;
  return laplace_base_rng(diff_likelihood<poisson_log_exposure_likelihood>(
                              L, y_and_ye, n_samples, msgs),
                          covariance_function, eta_dummy, theta_0, data_tuple,
                          rng, msgs, tolerance, max_num_steps,
                          hessian_block_size, solver,
                          max_steps_line_search, std::forward<Args>(args)...);
}
}  // namespace math
}  // namespace stan

#endif
