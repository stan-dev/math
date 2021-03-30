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
  template <typename T0, typename T1, typename K>
  T1 laplace_marginal_bernoulli_logit_lpmf
    (const std::vector<int>& y,
     const std::vector<int>& n_samples,
     const K& covariance_function,
     const Eigen::Matrix<T1, Eigen::Dynamic, 1>& phi,
     const std::vector<Eigen::VectorXd>& x,
     const std::vector<double>& delta,
     const std::vector<int>& delta_int,
     const Eigen::Matrix<T0, Eigen::Dynamic, 1>& theta_0,
     std::ostream* msgs = nullptr,
     double tolerance = 1e-6,
     long int max_num_steps = 100) {
    // TODO: change this to a VectorXd once we have operands & partials.
    Eigen::Matrix<T1, Eigen::Dynamic, 1> eta_dummy(0);
    return laplace_marginal_density(
      diff_bernoulli_logit(to_vector(n_samples), to_vector(y)),
      covariance_function,
      phi, eta_dummy, x, delta, delta_int,
      theta_0, msgs, tolerance, max_num_steps);
  }

  // Add signature that takes x as a matrix instead of a vector.
  template <typename T0, typename T1, typename K>
  T1 laplace_marginal_bernoulli_logit_lpmf
    (const std::vector<int>& y,
     const std::vector<int>& n_samples,
     const K& covariance_function,
     const Eigen::Matrix<T1, Eigen::Dynamic, 1>& phi,
     const Eigen::MatrixXd& x,
     const std::vector<double>& delta,
     const std::vector<int>& delta_int,
     const Eigen::Matrix<T0, Eigen::Dynamic, 1>& theta_0,
     std::ostream* msgs = nullptr,
     double tolerance = 1e-6,
     long int max_num_steps = 100) {
  // TODO: change this to a VectorXd once we have operands & partials.
  Eigen::Matrix<T1, Eigen::Dynamic, 1> eta_dummy(0);
  return laplace_marginal_density(
    diff_bernoulli_logit(to_vector(n_samples), to_vector(y)),
    covariance_function,
    phi, eta_dummy, x, delta, delta_int,
    theta_0, msgs, tolerance, max_num_steps);
  }

}  // namespace math
}  // namespace stan

#endif
