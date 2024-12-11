#ifndef STAN_MATH_MIX_PROB_LAPLACE_BERNOULLI_LOGIT_RNG_HPP
#define STAN_MATH_MIX_PROB_LAPLACE_BERNOULLI_LOGIT_RNG_HPP

#include <stan/math/mix/functor/laplace_base_rng.hpp>
#include <stan/math/mix/functor/laplace_likelihood.hpp>
#include <stan/math/mix/prob/laplace_marginal_bernoulli_logit_lpmf.hpp>

namespace stan {
namespace math {

/**
 * In a latent gaussian model,
 *
 *   theta ~ Normal(theta | 0, Sigma(phi))
 *   y ~ pi(y | theta)
 *
 * return a multivariate normal random variate sampled
 * from the gaussian approximation of p(theta | y, phi),
 * where the likelihood is a Bernoulli with logit link.
 * @tparam CovarFun Type of structure for covariance function.
 * @tparam ThetaMatrix Type for latent Gaussian variable.
 * @tparam RNG Type of rng number.
 * @tparam TrainTuple Type for observed/training covariates for covariance fn.
 * @tparam PredTuple Type for predictive covariates for covariance fn.
 * @tparam Args Type for variadic arguments for likelihood function.
 * @param y Vector Vector of total number of trials with a positive outcome.
 * @param n_samples Vector of number of trials.
 * @param theta_0 Initial guess for mode of Laplace approximation.
 * @param covariance_function Covariance function.
 * @param train_tuple Observed/training covariates for covariance function.
 * @param pred_tuple Predictive covariates for covariance function.
 * @param tolerance Tolerared change in objective function for Laplace approx.
 * @param max_num_steps Max number of iterations of Newton solver for Laplace
 *                      approx.
 * @param hessian_block_size Size of blocks for Hessian of log likelihood w.r.t
 *                           latent Gaussian variables.
 * @param solver Type of Newton solver. Each corresponds to a distinct choice
 *               of B matrix (i.e. application SWM formula):
 *               1. computes square-root of negative Hessian.
 *               2. computes square-root of covariance matrix.
 *               3. computes no square-root and uses LU decomposition.
 * @param max_steps_line_search Number of steps after which the algorithm gives
 *                              up on doing a linesearch. If 0, no linesearch.
 * @param rng Rng number.
 * @param msgs Streaming message for covariance functions.
 * @param args Arguments for log likelihood function.
 *
 */
template <typename CovarFun, typename ThetaMatrix, class RNG,
          typename TrainTuple, typename PredTuple, typename CovarArgs,
          require_eigen_t<ThetaMatrix>* = nullptr>
inline Eigen::VectorXd  // CHECK -- right return type
laplace_marginal_tol_bernoulli_logit_rng(
    const std::vector<int>& y, const std::vector<int>& n_samples,
    const ThetaMatrix& theta_0, CovarFun&& covariance_function,
    CovarArgs&& covar_args, TrainTuple&& train_tuple, PredTuple&& pred_tuple,
    const double tolerance, const int64_t max_num_steps,
    const int hessian_block_size, const int solver,
    const int max_steps_line_search, RNG& rng, std::ostream* msgs) {
  laplace_options ops{hessian_block_size, solver, max_steps_line_search,
                      tolerance, max_num_steps};
  Eigen::Matrix<double, 0, 0> eta_dummy;
  return laplace_base_rng(
      bernoulli_logit_likelihood{},
      std::forward_as_tuple(to_vector(y), n_samples), covariance_function,
      eta_dummy, theta_0, rng, msgs, ops, std::forward<TrainTuple>(train_tuple),
      std::forward<PredTuple>(pred_tuple), std::forward<CovarArgs>(covar_args));
}

/**
 * In a latent gaussian model,
 *
 *   theta ~ Normal(theta | 0, Sigma(phi))
 *   y ~ pi(y | theta)
 *
 * return a multivariate normal random variate sampled
 * from the gaussian approximation of p(theta | y, phi),
 * where the likelihood is a Bernoulli with logit link.
 * @tparam CovarFun Type of structure for covariance function.
 * @tparam ThetaMatrix Type for latent Gaussian variable.
 * @tparam RNG Type of rng number.
 * @tparam TrainTuple Type for observed/training covariates for covariance fn.
 * @tparam PredTuple Type for predictive covariates for covariance fn.
 * @tparam Args Type for variadic arguments for likelihood function.
 * @param y Vector Vector of total number of trials with a positive outcome.
 * @param n_samples Vector of number of trials.
 * @param theta_0 Initial guess for mode of Laplace approximation.
 * @param covariance_function Covariance function.
 * @param train_tuple Observed/training covariates for covariance function.
 * @param pred_tuple Predictive covariates for covariance function.
 * @param tolerance Tolerared change in objective function for Laplace approx.
 * @param max_num_steps Max number of iterations of Newton solver for Laplace
 *                      approx.
 * @param hessian_block_size Size of blocks for Hessian of log likelihood w.r.t
 *                           latent Gaussian variables.
 * @param solver Type of Newton solver. Each corresponds to a distinct choice
                 of B matrix (i.e. application SWM formula):
 *               1. computes square-root of negative Hessian.
 *               2. computes square-root of covariance matrix.
 *               3. computes no square-root and uses LU decomposition.
 * @param max_steps_line_search Number of steps after which the algorithm gives
 *                              up on doing a linesearch. If 0, no linesearch.
 * @param rng Rng number.
 * @param msgs Streaming message for covariance and likelihood functions.
 * @param args Arguments for log likelihood function.
 *
 */
template <typename CovarFun, typename ThetaMatrix, class RNG,
          typename TrainTuple, typename PredTuple, typename CovarArgs,
          require_eigen_t<ThetaMatrix>* = nullptr>
inline Eigen::VectorXd  // CHECK -- right return type
laplace_marginal_bernoulli_logit_rng(
    const std::vector<int>& y, const std::vector<int>& n_samples,
    const ThetaMatrix& theta_0, CovarFun&& covariance_function,
    CovarArgs&& covar_args, TrainTuple&& train_tuple, PredTuple&& pred_tuple,
    RNG& rng, std::ostream* msgs) {
  constexpr laplace_options ops{1, 1, 0, 1e-6, 100};
  Eigen::Matrix<double, 0, 0> eta_dummy;
  return laplace_base_rng(bernoulli_logit_likelihood{},
                          std::forward_as_tuple(to_vector(y), n_samples),
                          covariance_function, eta_dummy, theta_0, ops,
                          std::make_tuple(), std::make_tuple(), rng, msgs,
                          std::forward<CovarArgs>(covar_args));
}

}  // namespace math
}  // namespace stan

#endif
