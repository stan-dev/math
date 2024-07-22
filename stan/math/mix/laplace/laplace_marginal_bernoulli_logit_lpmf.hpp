#ifndef STAN_MATH_LAPLACE_LAPLACE_MARGINAL_BERNOULLI_LPMF_HPP
#define STAN_MATH_LAPLACE_LAPLACE_MARGINAL_BERNOULLI_LPMF_HPP

#include <stan/math/mix/laplace/laplace_marginal.hpp>
#include <stan/math/mix/laplace/laplace_likelihood_bernoulli_logit.hpp>

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
template <typename CovarF, typename ThetaMatrix, typename... Args,
          require_eigen_t<ThetaMatrix>* = nullptr>
inline auto laplace_marginal_tol_bernoulli_logit_lpmf(
    const std::vector<int>& y, const std::vector<int>& n_samples,
    double tolerance, long int max_num_steps, const int hessian_block_size,
    const int solver, const int max_steps_line_search,
    const ThetaMatrix& theta_0, CovarF&& covariance_function,
    std::ostream* msgs, Args&&... args) {
  // TODO: change this to a VectorXd once we have operands & partials.
  Eigen::Matrix<double, 0, 0> eta_dummy;
  return laplace_marginal_density(
      diff_likelihood<bernoulli_logit_likelihood>(
          bernoulli_logit_likelihood{}, to_vector(y), n_samples, msgs),
      covariance_function, eta_dummy, theta_0, msgs, tolerance, max_num_steps,
      hessian_block_size, solver, max_steps_line_search,
      std::forward<Args>(args)...);
}

template <typename CovarF, typename ThetaMatrix, typename... Args,
          require_eigen_t<ThetaMatrix>* = nullptr>
inline auto laplace_marginal_bernoulli_logit_lpmf(
    const std::vector<int>& y, const std::vector<int>& n_samples,
    const ThetaMatrix& theta_0, CovarF&& covariance_function,
    std::ostream* msgs, Args&&... args) {
  constexpr double tolerance = 1e-6;
  constexpr long int max_num_steps = 100;
  constexpr int hessian_block_size = 1;
  constexpr int solver = 1;
  constexpr int max_steps_line_search = 0;
  // TODO: change this to a VectorXd once we have operands & partials.
  Eigen::Matrix<double, 0, 0> eta_dummy;
  return laplace_marginal_density(
      diff_likelihood<bernoulli_logit_likelihood>(
          bernoulli_logit_likelihood{}, to_vector(y), n_samples, msgs),
      covariance_function, eta_dummy, theta_0, msgs, tolerance, max_num_steps,
      hessian_block_size, solver, max_steps_line_search,
      std::forward<Args>(args)...);
}

}  // namespace math
}  // namespace stan

#endif
