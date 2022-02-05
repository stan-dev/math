#ifndef STAN_MATH_LAPLACE_LAPLACE_MARGINAL_BERNOULLI_LPMF_HPP
#define STAN_MATH_LAPLACE_LAPLACE_MARGINAL_BERNOULLI_LPMF_HPP

#include <stan/math/laplace/laplace_marginal.hpp>
#include <stan/math/laplace/laplace_likelihood_bernoulli_logit.hpp>

namespace stan {
namespace math {
// EXPERIMENT
// Use the squared exponential kernel, for the time defined
// in the laplace_likelihood folder.
// In the final version, the user will provide the covariance
// function.

/**
 * Wrapper function around the laplace_marginal function for
 * a logistic Bernoulli likelihood. Returns the marginal density
 * p(y | phi) by marginalizing out the latent gaussian variable,
 * with a Laplace approximation. See the laplace_marginal function
 * for more details.
 *
 * @tparam T0 The type of the initial guess, theta_0.
 * @tparam T1 The type for the global parameter, phi.
 * @param[in] theta_0 the initial guess for the Laplace approximation.
 * @param[in] phi model parameters for the covariance function.
 * @param[in] x data for the covariance function.
 * @param[in] n_samples number of samples per group. First sufficient
 *            statistics.
 * @param[in] y total counts per group. Second sufficient statistics.
 * @param[in] tolerance controls the convergence criterion when finding
 *            the mode in the Laplace approximation.
 * @param[in] max_num_steps maximum number of steps before the Newton solver
 *            breaks and returns an error.
 */
template <typename T0, typename CovarF, typename... Args>
inline auto laplace_marginal_bernoulli_logit_lpmf(
    const std::vector<int>& y, const std::vector<int>& n_samples,
    CovarF&& covariance_function,
    const Eigen::Matrix<T0, Eigen::Dynamic, 1>& theta_0,
    std::ostream* msgs = nullptr, double tolerance = 1e-6,
    long int max_num_steps = 100, const int hessian_block_size = 0,
    const int solver = 1, const int do_line_search = 0,
    const int max_steps_line_search = 10, Args&&... args) {
  // TODO: change this to a VectorXd once we have operands & partials.
  Eigen::Matrix<double, Eigen::Dynamic, 1> eta_dummy(0);
  return laplace_marginal_density(
      diff_likelihood<bernoulli_logit_likelihood>(
          bernoulli_logit_likelihood{}, to_vector(y), n_samples, msgs),
      covariance_function, eta_dummy, theta_0, msgs, tolerance, max_num_steps,
      hessian_block_size, solver, do_line_search, max_steps_line_search,
      std::forward<Args>(args)...);
}
/*
// Add signature that takes x as a matrix instead of a vector.
template <typename T0, typename CovarF, typename... Args>
inline auto laplace_marginal_bernoulli_logit_lpmf(
    const std::vector<int>& y, const std::vector<int>& n_samples,
    CovarF&& covariance_function,
    const Eigen::Matrix<T0, Eigen::Dynamic, 1>& theta_0,
    std::ostream* msgs = nullptr, double tolerance = 1e-6,
    long int max_num_steps = 100, const int hessian_block_size = 0,
    const int solver = 1, const int do_line_search = 0,
    const int max_steps_line_search = 10, Args&&... args) {
  // TODO: change this to a VectorXd once we have operands & partials.
  Eigen::Matrix<double, Eigen::Dynamic, 1> eta_dummy(0);
  return laplace_marginal_density(
      diff_likelihood<bernoulli_logit_likelihood>(
          bernoulli_logit_likelihood{}, to_vector(y), n_samples, msgs),
      covariance_function, eta_dummy, theta_0, msgs, tolerance, max_num_steps,
     hessian_block_size, solver, do_line_search, max_steps_line_search,
std::forward<Args>(args)...);
}
*/
}  // namespace math
}  // namespace stan

#endif
