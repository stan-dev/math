#ifndef STAN_MATH_MIX_PROB_LAPLACE_MARGINAL_NEG_BINOMIAL_2_LOG_LPMF_HPP
#define STAN_MATH_MIX_PROB_LAPLACE_MARGINAL_NEG_BINOMIAL_2_LOG_LPMF_HPP

#include <stan/math/mix/functor/laplace_likelihood.hpp>
#include <stan/math/mix/functor/laplace_marginal_density.hpp>

namespace stan {
namespace math {

struct neg_binomial_2_log_likelihood {
  template <typename T_theta, typename T_eta, typename Y_t, typename Sums_t,
            typename NSamples>
  inline return_type_t<T_theta, T_eta> operator()(
      const Eigen::Matrix<T_theta, Eigen::Dynamic, 1>& theta,
      const Eigen::Matrix<T_eta, Eigen::Dynamic, 1>& eta, Y_t&& y,
      Sums_t&& sums, NSamples&& n_samples) const {
    T_eta eta_scalar = eta(0);
    return_type_t<T_theta, T_eta> logp = 0;
    for (size_t i = 0; i < y.size(); i++) {
      logp += binomial_coefficient_log(y(i) + eta_scalar - 1, y(i));
    }
    // CHECK -- is it better to vectorize this loop?
    Eigen::Matrix<T_theta, Eigen::Dynamic, 1> exp_theta = exp(theta);
    for (Eigen::Index i = 0; i < theta.size(); i++) {
      return_type_t<T_theta, T_eta> log_eta_plus_exp_theta
          = log(eta_scalar + exp_theta(i));
      logp += sums(i) * (theta(i) - log_eta_plus_exp_theta)
              + n_samples(i) * eta_scalar
                    * (log(eta_scalar) - log_eta_plus_exp_theta);
    }
    return logp;
  }
};

/**
 * Wrapper function around the laplace_marginal function for
 * a negative binomial likelihood. Uses the 2nd parameterization.
 * Returns the marginal density p(y | phi) by marginalizing
 * out the latent gaussian variable, with a Laplace approximation.
 * See the laplace_marginal function for more details.
 *
 * @tparam CovarFun The type of the initial guess, theta_0.
 * @tparam Eta The type for the global parameter, phi.
 * @tparam Theta0 The type of the initial guess, theta_0.
 * @tparam Args The type of arguments for the covariance function.
 * @param[in] y observations.
 * @param[in] y_index group to which each observation belongs. Each group
 *            is parameterized by one element of theta.
 * @param[in] n_samples
 * @param[in] sums
 * @param[in] eta non-marginalized model parameters for the likelihood.
 * @param[in] theta_0 the initial guess for the Laplace approximation.
 * @param[in] covariance_function a function which returns the prior covariance.
 * @param[in] tolerance controls the convergence criterion when finding
 *            the mode in the Laplace approximation.
 * @param[in] max_num_steps maximum number of steps before the Newton solver
 *            breaks and returns an error.
 * @param[in] hessian_block_size the size of the block for a block-diagonal
 *              Hessian of the log likelihood. If 0, the Hessian is stored
 *              inside a vector. If the Hessian is dense, this should be the
 *              size of the Hessian.
 * @param[in] solver Type of Newton solver. Each corresponds to a distinct
 *               choice of B matrix (i.e. application SWM formula):
 *               1. computes square-root of negative Hessian.
 *               2. computes square-root of covariance matrix.
 *               3. computes no square-root and uses LU decomposition.
 * @param[in] max_steps_line_search Number of steps after which the algorithm
 *                          gives up on doing a linesearch. If 0, no linesearch.
 * @param[in] msgs message stream for the covariance and likelihood function.
 * @param[in] args model parameters and data for the covariance functor.
template <typename CovarFun, typename Eta, typename Theta0, typename CovarArgs>
inline auto laplace_marginal_tol_neg_binomial_2_log_lpmf(
    const std::vector<int>& y, const std::vector<int>& y_index,
    const Eigen::VectorXd& n_samples, const Eigen::VectorXd& sums,
    const Eta& eta, const Theta0& theta_0, CovarFun&& covariance_function,
    CovarArgs&& covar_args, double tolerance, int64_t max_num_steps,
    const int hessian_block_size, const int solver,
    const int max_steps_line_search, std::ostream* msgs) {
  laplace_options ops{hessian_block_size, solver, max_steps_line_search,
                      tolerance, max_num_steps};
  return laplace_marginal_density(
      neg_binomial_2_log_likelihood{},
      std::forward_as_tuple(to_vector(y), y_index, n_samples, sums), eta,
      theta_0, std::forward<CovarFun>(covariance_function),
      std::forward<CovarArgs>(covar_args), ops, msgs);
}
 */

/**
 * Wrapper function around the laplace_marginal function for
 * a negative binomial likelihood. Uses the 2nd parameterization.
 * Returns the marginal density p(y | phi) by marginalizing
 * out the latent gaussian variable, with a Laplace approximation.
 * See the laplace_marginal function for more details.
 *
 * @tparam CovarFun The type of the initial guess, theta_0.
 * @tparam Theta0 The type of the initial guess, theta_0.
 * @tparam Eta The type of parameter arguments for the likelihood function.
 * @tparam Args Type of variadic arguments for covariance function.
 * @param[in] y observations.
 * @param[in] y_index group to which each observation belongs. Each group
 *            is parameterized by one element of theta.
 * @param[in] n_samples Number of count observations per group.
 * @param[in] sums Total number of counts per group.
 * @param[in] eta Parameter argument for likelihood function.
 * @param[in] theta_0 the initial guess for the Laplace approximation.
 * @param[in] covariance_function a function which returns the prior covariance.
 * @param[in] msgs  message stream for the covariance and likelihood function.
 * @param[in] args model parameters and data for the covariance functor.
template <typename CovarFun, typename Eta, typename Theta0, typename CovarArgs>
inline auto laplace_marginal_neg_binomial_2_log_lpmf(
    const std::vector<int>& y, const std::vector<int>& y_index,
    const Eigen::VectorXd& n_samples, const Eigen::VectorXd& sums,
    const Eta& eta, const Theta0& theta_0, CovarFun&& covariance_function,
    CovarArgs&& covar_args, std::ostream* msgs) {
  constexpr laplace_options ops{1, 1, 0, 1e-6, 100};
  return laplace_marginal_density(
      neg_binomial_2_log_likelihood{},
      std::forward_as_tuple(to_vector(y), y_index, n_samples, sums), eta,
      theta_0, std::forward<CovarFun>(covariance_function),
      std::forward<CovarArgs>(covar_args), ops, msgs);
}
 */
}  // namespace math
}  // namespace stan

#endif
