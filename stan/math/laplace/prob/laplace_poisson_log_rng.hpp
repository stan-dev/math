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
template <typename CovarFun, typename ThetaMatrix, class RNG,
          typename TrainTuple, typename PredTuple, typename... Args,
          require_eigen_t<ThetaMatrix>* = nullptr>
inline Eigen::VectorXd laplace_marginal_tol_poisson_log_rng(
    const std::vector<int>& y, const std::vector<int>& n_samples,
    const double tolerance, const long int max_num_steps,
    const int hessian_block_size, const int solver,
    const int max_steps_line_search, const ThetaMatrix& theta_0,
    CovarFun&& covariance_function, RNG& rng, std::ostream* msgs,
    TrainTuple&& train_tuple, PredTuple&& pred_tuple, Args&&... args) {
  Eigen::VectorXd eta_dummy;
  poisson_log_likelihood L;
  return laplace_base_rng(
      diff_likelihood<poisson_log_likelihood>(L, to_vector(y), n_samples, msgs),
      covariance_function, eta_dummy, theta_0, rng, msgs, tolerance,
      max_num_steps, hessian_block_size, solver, max_steps_line_search,
      std::forward<TrainTuple>(train_tuple),
      std::forward<PredTuple>(pred_tuple), std::forward<Args>(args)...);
}

/**
 * Overload for case where user passes exposure, ye.
 */
template <typename CovarFun, typename ThetaMatrix, class RNG,
          typename TrainTuple, typename PredTuple, typename... Args,
          require_eigen_t<ThetaMatrix>* = nullptr>
inline Eigen::VectorXd  // CHECK -- right return type
laplace_marginal_tol_poisson_2_log_rng(
    const std::vector<int>& y, const std::vector<int>& n_samples,
    const Eigen::VectorXd& ye, const double tolerance,
    const long int max_num_steps, const int hessian_block_size,
    const int solver, const int max_steps_line_search,
    const ThetaMatrix& theta_0, CovarFun&& covariance_function, RNG& rng,
    std::ostream* msgs, TrainTuple&& train_tuple, PredTuple&& pred_tuple,
    Args&&... args) {
  Eigen::VectorXd eta_dummy;
  Eigen::VectorXd y_vec = to_vector(y);
  Eigen::VectorXd y_and_ye(y_vec.size() + ye.size());
  y_and_ye << y_vec, ye;
  poisson_log_exposure_likelihood L;
  return laplace_base_rng(
      diff_likelihood<poisson_log_exposure_likelihood>(L, y_and_ye, n_samples,
                                                       msgs),
      covariance_function, eta_dummy, theta_0, rng, msgs, tolerance,
      max_num_steps, hessian_block_size, solver, max_steps_line_search,
      std::forward<TrainTuple>(train_tuple),
      std::forward<PredTuple>(pred_tuple), std::forward<Args>(args)...);
}

template <typename CovarFun, typename ThetaMatrix, class RNG,
          typename TrainTuple, typename PredTuple, typename... Args,
          require_eigen_t<ThetaMatrix>* = nullptr>
inline Eigen::VectorXd laplace_marginal_poisson_log_rng(
    const std::vector<int>& y, const std::vector<int>& n_samples,
    const ThetaMatrix& theta_0, CovarFun&& covariance_function, RNG& rng,
    std::ostream* msgs, TrainTuple&& train_tuple, PredTuple&& pred_tuple,
    Args&&... args) {
  Eigen::VectorXd eta_dummy;
  constexpr double tolerance = 1e-6;
  constexpr long int max_num_steps = 100;
  constexpr int hessian_block_size = 1;
  constexpr int solver = 1;
  constexpr int max_steps_line_search = 0;
  return laplace_base_rng(
      diff_likelihood<poisson_log_likelihood>(poisson_log_likelihood{},
                                              to_vector(y), n_samples, msgs),
      covariance_function, eta_dummy, theta_0, rng, msgs, tolerance,
      max_num_steps, hessian_block_size, solver, max_steps_line_search,
      std::forward<TrainTuple>(train_tuple),
      std::forward<PredTuple>(pred_tuple), std::forward<Args>(args)...);
}

/**
 * Overload for case where user passes exposure, ye.
 */
template <typename CovarFun, typename ThetaMatrix, class RNG,
          typename TrainTuple, typename PredTuple, typename... Args,
          require_eigen_t<ThetaMatrix>* = nullptr>
inline Eigen::VectorXd  // CHECK -- right return type
laplace_marginal_poisson_2_log_rng(const std::vector<int>& y,
                                   const std::vector<int>& n_samples,
                                   const Eigen::VectorXd& ye,
                                   const ThetaMatrix& theta_0,
                                   CovarFun&& covariance_function, RNG& rng,
                                   std::ostream* msgs, TrainTuple&& train_tuple,
                                   PredTuple&& pred_tuple, Args&&... args) {
  constexpr double tolerance = 1e-6;
  constexpr long int max_num_steps = 100;
  constexpr int hessian_block_size = 1;
  constexpr int solver = 1;
  constexpr int max_steps_line_search = 0;
  Eigen::VectorXd eta_dummy;
  Eigen::VectorXd y_vec = to_vector(y);
  Eigen::VectorXd y_and_ye(y_vec.size() + ye.size());
  y_and_ye << y_vec, ye;
  return laplace_base_rng(
      diff_likelihood<poisson_log_exposure_likelihood>(
          poisson_log_exposure_likelihood{}, y_and_ye, n_samples, msgs),
      covariance_function, eta_dummy, theta_0, rng, msgs, tolerance,
      max_num_steps, hessian_block_size, solver, max_steps_line_search,
      std::forward<TrainTuple>(train_tuple),
      std::forward<PredTuple>(pred_tuple), std::forward<Args>(args)...);
}

}  // namespace math
}  // namespace stan

#endif
