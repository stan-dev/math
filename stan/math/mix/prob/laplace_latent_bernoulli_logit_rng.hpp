#ifndef STAN_MATH_MIX_PROB_LAPLACE_LATENT_BERNOULLI_LOGIT_RNG_HPP
#define STAN_MATH_MIX_PROB_LAPLACE_LATENT_BERNOULLI_LOGIT_RNG_HPP

#include <stan/math/mix/functor/laplace_base_rng.hpp>
#include <stan/math/mix/functor/laplace_likelihood.hpp>
#include <stan/math/mix/prob/laplace_marginal_bernoulli_logit_lpmf.hpp>

namespace stan {
namespace math {

/**
 * In a latent gaussian model,
 *
 *   theta ~ Normal(theta|0, Sigma(phi))
 *   y ~ pi(y|theta)
 *
 * return a multivariate normal random variate sampled
 * from the gaussian approximation of p(theta | y, phi),
 * where the likelihood is a Bernoulli with logit link.
 * @tparam ThetaVec A type inheriting from `Eigen::EigenBase` with dynamic
 * sized rows and 1 column.
 * @tparam CovarFun A functor with an
 *  `operator()(CovarArgsElements...)` method. The `operator()` method should
 * accept as arguments the inner elements of `CovarArgs`. The return type of the
 * `operator()` method should be a type inheriting from `Eigen::EigenBase` with
 * dynamic sized rows and columns.
 * @tparam CovarArgs A tuple of types to passed as the first arguments of
 * `CovarFun::operator()`
 * @tparam RNG A valid boost rng type
 * @param y Vector Vector of total number of trials with a positive outcome.
 * @param n_samples Vector of number of trials.
 * @param theta_0 Initial guess for mode of Laplace approximation.
 * @param covariance_function Covariance function.
 * @param covar_args Observed/training covariates for covariance function.
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
 */
template <typename ThetaVec, typename CovarFun, typename CovarArgs,
          typename RNG, require_eigen_t<ThetaVec>* = nullptr>
inline Eigen::VectorXd  // CHECK -- right return type
laplace_latent_tol_bernoulli_logit_rng(
    const std::vector<int>& y, const std::vector<int>& n_samples,
    ThetaVec&& theta_0, CovarFun&& covariance_function, CovarArgs&& covar_args,
    const double tolerance, const int64_t max_num_steps,
    const int hessian_block_size, const int solver,
    const int max_steps_line_search, RNG& rng, std::ostream* msgs) {
  laplace_options ops{hessian_block_size, solver, max_steps_line_search,
                      tolerance, max_num_steps};
  return laplace_base_rng(bernoulli_logit_likelihood{},
                          std::forward_as_tuple(to_vector(y), n_samples),
                          std::forward<ThetaVec>(theta_0),
                          std::forward<CovarFun>(covariance_function),
                          std::forward<CovarArgs>(covar_args), ops, rng, msgs);
}

/**
 * In a latent gaussian model,
 *
 *   theta ~ Normal(theta|0, Sigma(phi))
 *   y ~ pi(y|theta)
 *
 * return a multivariate normal random variate sampled
 * from the gaussian approximation of p(theta | y, phi),
 * where the likelihood is a Bernoulli with logit link.
 * @tparam ThetaVec A type inheriting from `Eigen::EigenBase` with dynamic
 * sized rows and 1 column.
 * @tparam CovarFun A functor with an
 *  `operator()(CovarArgsElements..., {TrainTupleElements...|
 PredTupleElements...})`
 *  method. The `operator()` method should accept as arguments the
 *  inner elements of `CovarArgs`. The return type of the `operator()` method
 *  should be a type inheriting from `Eigen::EigenBase` with dynamic sized
 *  rows and columns.
 * @tparam CovarArgs A tuple of types to passed as the first arguments of
 `CovarFun::operator()`
 * @tparam RNG A valid boost rng type
 * @param y Vector Vector of total number of trials with a positive outcome.
 * @param n_samples Vector of number of trials.
 * @param theta_0 Initial guess for mode of Laplace approximation.
 * @param covariance_function Covariance function.
 * @param covar_args Observed/training covariates for covariance function.
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
 */
template <typename CovarFun, typename ThetaVec, typename CovarArgs,
          typename RNG,
          require_eigen_t<ThetaVec>* = nullptr>
inline Eigen::VectorXd  // CHECK -- right return type
laplace_latent_bernoulli_logit_rng(const std::vector<int>& y,
                                   const std::vector<int>& n_samples,
                                   ThetaVec&& theta_0,
                                   CovarFun&& covariance_function,
                                   CovarArgs&& covar_args, RNG& rng,
                                   std::ostream* msgs) {
  constexpr laplace_options ops{1, 1, 0, 1e-6, 100};
  return laplace_base_rng(bernoulli_logit_likelihood{},
                          std::forward_as_tuple(to_vector(y), n_samples),
                          std::forward<ThetaVec>(theta_0),
                          std::forward<CovarFun>(covariance_function),
                          std::forward<CovarArgs>(covar_args), ops, rng, msgs);
}

}  // namespace math
}  // namespace stan

#endif
