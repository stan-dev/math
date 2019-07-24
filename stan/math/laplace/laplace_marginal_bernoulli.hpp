#ifndef STAN_MATH_LAPLACE_LAPLACE_MARGINAL_BERNOULLI_HPP
#define STAN_MATH_LAPLACE_LAPLACE_MARGINAL_BERNOULLI_HPP

#include <stan/math/laplace/laplace_marginal.hpp>
#include <stan/math/laplace/laplace_likelihood.hpp>

namespace stan {
namespace math {
  // Use the squared exponential curvature, for the time defined
  // in the laplace_likelihood folder.
  // In the final version, the user will provide the covariance
  // function.

  /**
   * Wrapper function to be exposed to the Stan language.
   * This first prototype will enforce the squared kernel for the
   * covariance matrix. Future versions will require the user to
   * pass their covariance function.
   */
  template <typename T0, typename T1>
  T1 laplace_marginal_bernoulli 
               (const Eigen::Matrix<T0, Eigen::Dynamic, 1>& theta_0,
                const Eigen::Matrix<T1, Eigen::Dynamic, 1>& phi,
                const std::vector<Eigen::VectorXd>& x,
                const std::vector<int>& n_samples,
                const std::vector<int>& y,
                // const K& covariance_function,
                double tolerance = 1e-6,
                long int max_num_steps = 100) {
    return laplace_marginal_density(theta_0, phi, x,
       diff_logistic_log(to_vector(n_samples), to_vector(y)),
       sqr_exp_kernel_functor(),
       tolerance, max_num_steps);
  }

}  // namespace math
}  // namespace stan

#endif
