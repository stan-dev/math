#ifndef STAN_MATH_MIX_PROB_LAPLACE_POISSON_LOG_EXPOSURE_RNG_HPP
#define STAN_MATH_MIX_PROB_LAPLACE_POISSON_LOG_EXPOSURE_RNG_HPP

#include <stan/math/mix/functor/laplace_base_rng.hpp>
#include <stan/math/mix/functor/laplace_likelihood.hpp>
#include <stan/math/mix/prob/laplace_marginal_poisson_log_exposure_lpmf.hpp>

namespace stan {
namespace math {

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
  laplace_options ops{hessian_block_size, solver,
    max_steps_line_search, tolerance, max_num_steps};
  return laplace_base_rng(poisson_log_exposure_likelihood{},
      std::forward_as_tuple(y_and_ye, n_samples),
      covariance_function, eta_dummy, theta_0, rng, msgs, ops,
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
  constexpr laplace_options ops{1, 1, 0, 1e-6, 100};
  return laplace_base_rng(poisson_log_exposure_likelihood{},
      std::forward_as_tuple(y_and_ye, n_samples),
      covariance_function, eta_dummy, theta_0, rng, msgs, ops,
      std::forward<TrainTuple>(train_tuple),
      std::forward<PredTuple>(pred_tuple), std::forward<Args>(args)...);
}

}  // namespace math
}  // namespace stan

#endif
