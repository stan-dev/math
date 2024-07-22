#ifndef STAN_MATH_LAPLACE_LAPLACE_MARGINAL_POISSON_LOG_LPMF_HPP
#define STAN_MATH_LAPLACE_LAPLACE_MARGINAL_POISSON_LOG_LPMF_HPP

#include <stan/math/mix/laplace/laplace_marginal.hpp>
#include <stan/math/mix/laplace/laplace_likelihood_poisson_log.hpp>
#include <stan/math/mix/laplace/laplace_likelihood_general.hpp>

namespace stan {
namespace math {
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
    double tolerance, long int max_num_steps, const int hessian_block_size,
    const int solver, const int max_steps_line_search, const ThetaVec& theta_0,
    CovarFun&& covariance_function, std::ostream* msgs, Args&&... args) {
  // TODO: change this to a VectorXd once we have operands & partials.
  Eigen::Matrix<double, 0, 0> eta_dummy;
  return laplace_marginal_density(
      diff_likelihood<poisson_log_likelihood>(poisson_log_likelihood{},
                                              to_vector(y), n_samples, msgs),
      covariance_function, eta_dummy, theta_0, msgs, tolerance, max_num_steps,
      hessian_block_size, solver, max_steps_line_search, args...);
}

template <typename CovarFun, typename YeVec, typename ThetaVec,
          typename... Args,
          require_all_eigen_vector_t<YeVec, ThetaVec>* = nullptr>
inline auto laplace_marginal_tol_poisson_2_log_lpmf(
    const std::vector<int>& y, const std::vector<int>& n_samples,
    const YeVec& ye, double tolerance, long int max_num_steps,
    const int hessian_block_size, const int solver,
    const int max_steps_line_search, const ThetaVec& theta_0,
    CovarFun&& covariance_function, std::ostream* msgs, Args&&... args) {
  // TODO: change this to a VectorXd once we have operands & partials.
  Eigen::Matrix<double, 0, 0> eta_dummy;
  Eigen::VectorXd y_vec = to_vector(y);
  Eigen::VectorXd y_and_ye(y_vec.size() + ye.size());
  y_and_ye << y_vec, ye;
  return laplace_marginal_density(
      diff_likelihood<poisson_log_exposure_likelihood>(
          poisson_log_exposure_likelihood{}, y_and_ye, n_samples, msgs),
      std::forward<CovarFun>(covariance_function), eta_dummy, theta_0, msgs,
      tolerance, max_num_steps, hessian_block_size, solver,
      max_steps_line_search, std::forward<Args>(args)...);
}

template <typename CovarFun, typename ThetaVec, typename... Args,
          require_eigen_vector_t<ThetaVec>* = nullptr>
inline auto laplace_marginal_poisson_log_lpmf(const std::vector<int>& y,
                                              const std::vector<int>& n_samples,
                                              const ThetaVec& theta_0,
                                              CovarFun&& covariance_function,
                                              std::ostream* msgs,
                                              Args&&... args) {
  // TODO: change this to a VectorXd once we have operands & partials.
  Eigen::Matrix<double, 0, 0> eta_dummy;
  constexpr double tolerance = 1e-6;
  constexpr long int max_num_steps = 100;
  constexpr int hessian_block_size = 1;
  constexpr int solver = 1;
  constexpr int max_steps_line_search = 0;
  return laplace_marginal_density(
      diff_likelihood<poisson_log_likelihood>(poisson_log_likelihood{},
                                              to_vector(y), n_samples, msgs),
      covariance_function, eta_dummy, theta_0, msgs, tolerance, max_num_steps,
      hessian_block_size, solver, max_steps_line_search, args...);
}

template <typename CovarFun, typename YeVec, typename ThetaVec,
          typename... Args,
          require_all_eigen_vector_t<YeVec, ThetaVec>* = nullptr>
inline auto laplace_marginal_poisson_2_log_lpmf(
    const std::vector<int>& y, const std::vector<int>& n_samples,
    const YeVec& ye, const ThetaVec& theta_0, CovarFun&& covariance_function,
    std::ostream* msgs, Args&&... args) {
  // TODO: change this to a VectorXd once we have operands & partials.
  Eigen::Matrix<double, 0, 0> eta_dummy;
  Eigen::VectorXd y_vec = to_vector(y);
  Eigen::VectorXd y_and_ye(y_vec.size() + ye.size());
  y_and_ye << y_vec, ye;
  constexpr double tolerance = 1e-6;
  constexpr long int max_num_steps = 100;
  constexpr int hessian_block_size = 1;
  constexpr int solver = 1;
  constexpr int max_steps_line_search = 0;
  return laplace_marginal_density(
      diff_likelihood<poisson_log_exposure_likelihood>(
          poisson_log_exposure_likelihood{}, y_and_ye, n_samples, msgs),
      std::forward<CovarFun>(covariance_function), eta_dummy, theta_0, msgs,
      tolerance, max_num_steps, hessian_block_size, solver,
      max_steps_line_search, std::forward<Args>(args)...);
}

}  // namespace math
}  // namespace stan

#endif
