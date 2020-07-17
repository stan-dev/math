 #ifndef STAN_MATH_LAPLACE_LAPLACE_MARGINAL_BERNOULLI_HPP
#define STAN_MATH_LAPLACE_LAPLACE_MARGINAL_BERNOULLI_HPP

#include <stan/math/laplace/laplace_marginal.hpp>
#include <stan/math/laplace/laplace_likelihood.hpp>

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
  // TODO: deprecate the below function. No default functor.
  template <typename T0, typename T1>
  T1 laplace_marginal_bernoulli_logit
               (const std::vector<int>& y,
                const std::vector<int>& n_samples,
                // const K& covariance function,
                const Eigen::Matrix<T1, Eigen::Dynamic, 1>& phi,
                const std::vector<Eigen::VectorXd>& x,
                const std::vector<double>& delta,
                const std::vector<int>& delta_int,
                const Eigen::Matrix<T0, Eigen::Dynamic, 1>& theta_0,
                std::ostream* msgs = nullptr,
                double tolerance = 1e-6,
                long int max_num_steps = 100) {
    return laplace_marginal_density(
      diff_logistic_log(to_vector(n_samples), to_vector(y)),
      sqr_exp_kernel_functor(),
      phi, x, delta, delta_int,
      theta_0, msgs, tolerance, max_num_steps);
  }

  // Add signature that takes in a Kernel functor specified by the user.
  template <typename T0, typename T1, typename K>
  T1 laplace_marginal_bernoulli_logit
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
    return laplace_marginal_density(
      diff_logistic_log(to_vector(n_samples), to_vector(y)),
      covariance_function,
      phi, x, delta, delta_int,
      theta_0, msgs, tolerance, max_num_steps);
  }

  // Add signature that takes x as a matrix instead of a vector.
  template <typename T0, typename T1, typename K>
  T1 laplace_marginal_bernoulli_logit
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
  return laplace_marginal_density(
    diff_logistic_log(to_vector(n_samples), to_vector(y)),
    covariance_function,
    phi, x, delta, delta_int,
    theta_0, msgs, tolerance, max_num_steps);
  }

}  // namespace math
}  // namespace stan

#endif
