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
inline auto  // CHECK -- right return type
laplace_marginal_tol_poisson_2_log_rng(
    const std::vector<int>& y, const std::vector<int>& n_samples,
    const Eigen::VectorXd& ye, const ThetaMatrix& theta_0,
    CovarFun&& covariance_function, RNG& rng, TrainTuple&& train_tuple,
    PredTuple&& pred_tuple, const double tolerance, const int64_t max_num_steps,
    const int hessian_block_size, const int solver,
    const int max_steps_line_search, std::ostream* msgs, Args&&... args) {
  Eigen::Matrix<double, 0, 0> eta_dummy;
  laplace_options ops{hessian_block_size, solver, max_steps_line_search,
                      tolerance, max_num_steps};
  return laplace_base_rng(poisson_log_exposure_likelihood{},
                          std::forward_as_tuple(y, ye, n_samples),
                          covariance_function, eta_dummy, theta_0, ops,
                          std::forward<TrainTuple>(train_tuple),
                          std::forward<PredTuple>(pred_tuple), rng, msgs,
                          std::forward<Args>(args)...);
}

/**
 * Overload for case where user passes exposure, ye.
 */
template <typename CovarFun, typename ThetaMatrix, class RNG,
          typename TrainTuple, typename PredTuple, typename... Args,
          require_eigen_t<ThetaMatrix>* = nullptr>
inline auto  // TODO(Steve): Allow scalar or std vector return
laplace_marginal_poisson_2_log_rng(
    const std::vector<int>& y, const std::vector<int>& n_samples,
    const Eigen::VectorXd& ye, const ThetaMatrix& theta_0,
    CovarFun&& covariance_function, TrainTuple&& train_tuple,
    PredTuple&& pred_tuple, RNG& rng, std::ostream* msgs, Args&&... args) {
  Eigen::Matrix<double, 0, 0> eta_dummy;
  constexpr laplace_options ops{1, 1, 0, 1e-6, 100};
  return laplace_base_rng(poisson_log_exposure_likelihood{},
                          std::forward_as_tuple(to_vector(y), ye, n_samples),
                          covariance_function, eta_dummy, theta_0, ops,
                          std::forward<TrainTuple>(train_tuple),
                          std::forward<PredTuple>(pred_tuple), rng, msgs,
                          std::forward<Args>(args)...);
}

}  // namespace math
}  // namespace stan

#endif
