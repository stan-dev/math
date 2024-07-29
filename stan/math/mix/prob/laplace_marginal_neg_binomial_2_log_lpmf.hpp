#ifndef STAN_MATH_MIX_PROB_LAPLACE_MARGINAL_NEG_BINOMIAL_2_LOG_LPMF_HPP
#define STAN_MATH_MIX_PROB_LAPLACE_MARGINAL_NEG_BINOMIAL_2_LOG_LPMF_HPP

#include <stan/math/mix/functor/laplace_likelihood.hpp>
#include <stan/math/mix/functor/laplace_marginal_density.hpp>

namespace stan {
namespace math {

struct neg_binomial_2_log_likelihood {
  template <typename T_theta, typename T_eta, typename Y_t, typename Sums_t, typename NSamples>
  inline return_type_t<T_theta, T_eta> operator()(
      const Eigen::Matrix<T_theta, Eigen::Dynamic, 1>& theta,
      const Eigen::Matrix<T_eta, Eigen::Dynamic, 1>& eta,
      Y_t&& y, Sums_t&& sums, NSamples&& n_samples) const {
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
 * @tparam T0 The type of the initial guess, theta_0.
 * @tparam T1 The type for the global parameter, phi.
 * @param[in] y observations.
 * @param[in] y_index group to which each observation belongs. Each group
 *            is parameterized by one element of theta.
 * @param[in] covariance a function which returns the prior covariance.
 * @param[in] phi model parameters for the covariance functor.
 * @param[in] eta non-marginalized model parameters for the likelihood.
 * @param[in] x data for the covariance functor.
 * @param[in] delta additional real data for the covariance functor.
 * @param[in] delta_int additional integer data for covariance functor.
 * @param[in] theta_0 the initial guess for the Laplace approximation.
 * @param[in] tolerance controls the convergence criterion when finding
 *            the mode in the Laplace approximation.
 * @param[in] max_num_steps maximum number of steps before the Newton solver
 *            breaks and returns an error.
 */
template <typename CovarFun, typename Eta, typename Theta0, typename... Args>
inline auto laplace_marginal_tol_neg_binomial_2_log_lpmf(
    const std::vector<int>& y, const std::vector<int>& y_index,
    const Eigen::VectorXd& n_samples, const Eigen::VectorXd& sums,
    const Eta& eta, double tolerance, long int max_num_steps,
    const int hessian_block_size, const int solver,
    const int max_steps_line_search, const Theta0& theta_0,
    CovarFun&& covariance_function, std::ostream* msgs, Args&&... args) {
  laplace_options ops{hessian_block_size, solver, max_steps_line_search,
                      tolerance, max_num_steps};
  return laplace_marginal_density(neg_binomial_2_log_likelihood{},
      std::forward_as_tuple(to_vector(y), y_index, n_samples, sums),
      std::forward<CovarFun>(covariance_function), eta, theta_0, msgs,
      tolerance, max_num_steps, hessian_block_size, solver,
      max_steps_line_search, std::forward<Args>(args)...);
}

template <typename CovarFun, typename Eta, typename Theta0, typename... Args>
inline auto laplace_marginal_neg_binomial_2_log_lpmf(
    const std::vector<int>& y, const std::vector<int>& y_index,
    const Eigen::VectorXd& n_samples, const Eigen::VectorXd& sums,
    const Eta& eta, const Theta0& theta_0, CovarFun&& covariance_function,
    std::ostream* msgs, Args&&... args) {
  constexpr laplace_options ops{1, 1, 0, 1e-6, 100};
  return laplace_marginal_density(neg_binomial_2_log_likelihood{},
      std::forward_as_tuple(to_vector(y), y_index, n_samples, sums),
      std::forward<CovarFun>(covariance_function), eta, theta_0, msgs, ops,
      std::forward<Args>(args)...);
}
}  // namespace math
}  // namespace stan

#endif
