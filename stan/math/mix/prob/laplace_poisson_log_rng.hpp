#ifndef STAN_MATH_MIX_PROB_LAPLACE_POISSON_LOG_RNG_HPP
#define STAN_MATH_MIX_PROB_LAPLACE_POISSON_LOG_RNG_HPP

#include <stan/math/mix/functor/laplace_base_rng.hpp>
#include <stan/math/mix/functor/laplace_likelihood.hpp>
#include <stan/math/mix/prob/laplace_marginal_poisson_log_exposure_lpmf.hpp>
#include <stan/math/mix/prob/laplace_marginal_poisson_log_lpmf.hpp>

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
 * @tparam CovarFun
 * @tparam ThetaMatrix
 * @tparam RNG
 * @tparam TrainTuple
 * @tparam PredTuple
 * @tparam Args
 * @param y
 * @param n_samples
 * @param theta_0
 * @param covariance_function
 * @param train_tuple
 * @param pred_tuple
 * @param tolerance
 * @param max_num_steps
 * @param hessian_block_size
 * @param solver
 * @param max_steps_line_search
 * @param rng
 * @param msgs
 * @param args
 *
 */
template <typename CovarFun, typename ThetaMatrix, class RNG,
          typename TrainTuple, typename PredTuple, typename CovarArgs,
          require_eigen_t<ThetaMatrix>* = nullptr>
inline Eigen::VectorXd laplace_marginal_tol_poisson_log_rng(
    const std::vector<int>& y, const std::vector<int>& n_samples,
    const ThetaMatrix& theta_0, CovarFun&& covariance_function,
    CovarArgs&& covar_args,
    TrainTuple&& train_tuple, PredTuple&& pred_tuple, const double tolerance,
    const int64_t max_num_steps, const int hessian_block_size, const int solver,
    const int max_steps_line_search, RNG& rng, std::ostream* msgs) {
  Eigen::VectorXd eta_dummy;
  laplace_options ops{hessian_block_size, solver, max_steps_line_search,
                      tolerance, max_num_steps};
  return laplace_base_rng(poisson_log_likelihood{},
                          std::forward_as_tuple(to_vector(y), n_samples),
                          covariance_function, eta_dummy, theta_0, ops,
                          std::forward<TrainTuple>(train_tuple),
                          std::forward<PredTuple>(pred_tuple), rng, msgs,
                          std::forward<CovarArgs>(covar_args));
}

/**
 * In a latent gaussian model,
 *
 *   theta ~ Normal(theta | 0, Sigma(phi))
 *   y ~ pi(y | theta)
 *
 * return a multivariate normal random variate sampled
 * from the gaussian approximation of p(theta | y, phi)
 * where the likelihood is a Poisson with a log link.
 * @tparam CovarFun
 * @tparam ThetaMatrix
 * @tparam RNG
 * @tparam TrainTuple
 * @tparam PredTuple
 * @tparam Args
 * @param y
 * @param n_samples
 * @param theta_0
 * @param covariance_function
 * @param train_tuple
 * @param pred_tuple
 * @param rng
 * @param msgs
 * @param args
 *
 */
template <typename CovarFun, typename ThetaMatrix, class RNG,
          typename TrainTuple, typename PredTuple, typename CovarArgs,
          require_eigen_t<ThetaMatrix>* = nullptr>
inline Eigen::VectorXd laplace_marginal_poisson_log_rng(
    const std::vector<int>& y, const std::vector<int>& n_samples,
    const ThetaMatrix& theta_0, CovarFun&& covariance_function,
    CovarArgs&& covar_args,
    TrainTuple&& train_tuple, PredTuple&& pred_tuple, RNG& rng,
    std::ostream* msgs) {
  Eigen::VectorXd eta_dummy;
  constexpr laplace_options ops{1, 1, 0, 1e-6, 100};
  return laplace_base_rng(poisson_log_likelihood{},
                          std::forward_as_tuple(to_vector(y), n_samples),
                          covariance_function, eta_dummy, theta_0, ops,
                          std::forward<TrainTuple>(train_tuple),
                          std::forward<PredTuple>(pred_tuple), rng, msgs,
                          std::forward<CovarArgs>(covar_args));
}

}  // namespace math
}  // namespace stan

#endif
