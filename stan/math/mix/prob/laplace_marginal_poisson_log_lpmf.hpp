#ifndef STAN_MATH_MIX_PROB_LAPLACE_MARGINAL_POISSON_LOG_LPMF_HPP
#define STAN_MATH_MIX_PROB_LAPLACE_MARGINAL_POISSON_LOG_LPMF_HPP

#include <stan/math/mix/functor/laplace_likelihood.hpp>
#include <stan/math/mix/functor/laplace_marginal_density.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <Eigen/Sparse>

namespace stan {
namespace math {

struct poisson_log_likelihood {
  /**
   * Returns the lpmf for a Poisson with a log link across
   * multiple groups. No need to compute the log normalizing constant.
   * @tparam T_theta Type of the log Poisson rate.
   * @tparam T_eta Type of the auxiliary parameter (not used here).
   * @param[in] theta log Poisson rate for each group.
   * @param[in] y sum of counts in each group.
   * @param[in] delta_int number of observations in each group.
   * return lpmf for a Poisson with a log link.
   */
  template <typename Theta, typename Eta,
            require_eigen_vector_t<Theta>* = nullptr,
            require_eigen_t<Eta>* = nullptr>
  inline auto operator()(const Theta& theta, const Eta& /* eta */,
                         const Eigen::VectorXd& y,
                         const std::vector<int>& delta_int,
                         std::ostream* pstream) const {
    auto n_samples = to_vector(delta_int);
    return -sum(lgamma(add(y, 1))) + dot_product(theta, y)
           - dot_product(n_samples, exp(theta));
  }
};

/**
 * Wrapper function around the laplace_marginal function for
 * a log poisson likelihood. Returns the marginal density
 * p(y | phi) by marginalizing out the latent gaussian variable,
 * with a Laplace approximation. See the laplace_marginal function
 * for more details.
 *
 * @tparam T0 The type of the initial guess, theta_0.
 * @tparam T1 The type for the global parameter, phi.
 * @param[in] y total counts per group. Second sufficient statistics.
 * @param[in] n_samples number of samples per group. First sufficient
 *            statistics.
 * @param[in] covariance a function which returns the prior covariance.
 * @param[in] phi model parameters for the covariance functor.
 * @param[in] x data for the covariance functor.
 * @param[in] delta additional real data for the covariance functor.
 * @param[in] delta_int additional integer data for covariance functor.
 * @param[in] theta_0 the initial guess for the Laplace approximation.
 * @param[in] tolerance controls the convergence criterion when finding
 *            the mode in the Laplace approximation.
 * @param[in] max_num_steps maximum number of steps before the Newton solver
 *            breaks and returns an error.
 */
template <typename CovarFun, typename ThetaVec, typename... Args,
          require_all_eigen_vector_t<ThetaVec>* = nullptr>
inline auto laplace_marginal_tol_poisson_log_lpmf(
    const std::vector<int>& y, const std::vector<int>& n_samples,
    const ThetaVec& theta_0, CovarFun&& covariance_function, double tolerance,
    int64_t max_num_steps, const int hessian_block_size, const int solver,
    const int max_steps_line_search, std::ostream* msgs, Args&&... args) {
  Eigen::Matrix<double, 0, 0> eta_dummy;
  laplace_options ops{hessian_block_size, solver, max_steps_line_search,
                      tolerance, max_num_steps};
  return laplace_marginal_density(
      poisson_log_likelihood{}, std::forward_as_tuple(to_vector(y), n_samples),
      covariance_function, eta_dummy, theta_0, msgs, ops, args...);
}

template <typename CovarFun, typename ThetaVec, typename... Args,
          require_eigen_vector_t<ThetaVec>* = nullptr>
inline auto laplace_marginal_poisson_log_lpmf(const std::vector<int>& y,
                                              const std::vector<int>& n_samples,
                                              const ThetaVec& theta_0,
                                              CovarFun&& covariance_function,
                                              std::ostream* msgs,
                                              Args&&... args) {
  Eigen::Matrix<double, 0, 0> eta_dummy;
  constexpr laplace_options ops{1, 1, 0, 1e-6, 100};
  return laplace_marginal_density(
      poisson_log_likelihood{}, std::forward_as_tuple(to_vector(y), n_samples),
      covariance_function, eta_dummy, theta_0, msgs, ops, args...);
}

}  // namespace math
}  // namespace stan

#endif
