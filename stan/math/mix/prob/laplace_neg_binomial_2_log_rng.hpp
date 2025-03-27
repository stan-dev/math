#ifndef STAN_MATH_MIX_PROB_LAPLACE_NEG_BINOMIAL_2_LOG_RNG_HPP
#define STAN_MATH_MIX_PROB_LAPLACE_NEG_BINOMIAL_2_LOG_RNG_HPP

#include <stan/math/mix/functor/laplace_base_rng.hpp>
#include <stan/math/mix/functor/laplace_likelihood.hpp>

namespace stan {
namespace math {

/**
 * In a latent gaussian model,
 *
 *   theta ~ Normal(0, Sigma(phi))
 *   y ~ p(y|theta,phi)
 *
 * return a sample from the Laplace approximation to p(theta|y,phi).
 * The Laplace approximation is computed using a Newton solver.
 * In this specialized function, the likelihood p(y|theta) is a
 * Negative Binomial with a log link. This function uses the second
 * parameterization of the Negative Binomial.
 *
 * @tparam CovarFun A functor with an
 *  `operator()(CovarArgsElements..., {TrainTupleElements...|
 PredTupleElements...})`
 *  method. The `operator()` method should accept as arguments the
 *  inner elements of `CovarArgs`, followed by either the inner elements of
 *  `TrainTuple` or `PredTuple`. The return type of the `operator()` method
 *  should be a type inheriting from `Eigen::EigenBase` with dynamic sized
 *  rows and columns.
 * @tparam Eta A type for the overdispersion parameter.
 * @tparam ThetaMatrix A type inheriting from `Eigen::EigenBase` with dynamic
 * sized rows and 1 column.
 * @tparam RNG A valid boost rng type
 * @tparam TrainTuple A tuple of types to passed as the end arguments of
 `CovarFun::operator()`
 * @tparam PredTuple  A tuple of types to passed as the end arguments of
 `CovarFun::operator()`
 * @tparam CovarArgs A tuple of types to passed as the first arguments of
 `CovarFun::operator()`
 * @param y Observed counts.
 * @param y_index Index indicating which group each observation belongs to.
 * @param eta Overdisperison parameter.
 * @param theta_0 Initial guess for the Newton solver.
 * @param covariance_function Function that returns prior covariance matrix.
 * @param covar_args arguments for the covariance function.
 * @param train_tuple additional arguments for the covariance function,
 *                    e.g. covariates which correspond to observed data.
 * @param pred_tuple additional arguments for the covariance function,
 *                   e.g. covariates for out-of-sample data.
 * @param tolerance tolerated norm of the gradient when finding the mode of
                    p(theta|y,phi) for the Laplace approximation.
 * @param max_num_steps maximum number of steps before the Newton solver
 *                      breaks and returns an error.
 * @param hessian_block_size the size of the block for a block-diagonal
 *              Hessian of the log likelihood, i.e second derivative of
 *              log p(y|theta,phi) wrt theta.
 * @param solver Type of Newton solver. Each corresponds to a distinct choice
 *               of B matrix (i.e. application SWM formula):
 *               1. computes square-root of negative Hessian.
 *               2. computes square-root of covariance matrix.
 *               3. computes no square-root and uses LU decomposition.
 * @param max_steps_line_search Number of steps after which the algorithm
 *                        gives up on doing a linesearch. If 0, no linesearch.
 * @param rng seed for rng.
 * @param msgs message stream for the covariance and likelihood function.
 */
template <typename CovarFun, typename Eta, typename ThetaMatrix, class RNG,
          typename TrainTuple, typename PredTuple, typename CovarArgs,
          require_eigen_t<ThetaMatrix>* = nullptr>
inline Eigen::VectorXd laplace_marginal_tol_neg_binomial_2_log_rng(
    const std::vector<int>& y, const std::vector<int>& y_index, const Eta& eta,
    const ThetaMatrix& theta_0, CovarFun&& covariance_function,
    CovarArgs&& covar_args, TrainTuple&& train_tuple, PredTuple&& pred_tuple,
    const double tolerance, const int64_t max_num_steps,
    const int hessian_block_size, const int solver,
    const int max_steps_line_search, RNG& rng, std::ostream* msgs) {
  const laplace_options ops{hessian_block_size, solver, max_steps_line_search,
                      tolerance, max_num_steps};
  return laplace_base_rng(neg_binomial_2_log_likelihood{},
                          std::forward_as_tuple(eta, y, y_index), theta_0,
                          std::forward<CovarFun>(covariance_function),
                          std::forward<CovarArgs>(covar_args),
                          std::forward<TrainTuple>(train_tuple),
                          std::forward<PredTuple>(pred_tuple), ops, rng, msgs);
}

/**
 * In a latent gaussian model,
 *
 *   theta ~ Normal(0, Sigma(phi))
 *   y ~ p(y|theta,phi)
 *
 * return a sample from the Laplace approximation to p(theta|y,phi).
 * The Laplace approximation is computed using a Newton solver.
 * In this specialized function, the likelihood p(y|theta) is a
 * Negative Binomial with a log link. This function uses the second
 * parameterization of the Negative Binomial.
 *
 * @tparam CovarFun A functor with an
 *  `operator()(CovarArgsElements..., {TrainTupleElements...|
 * PredTupleElements...})` method. The `operator()` method should accept as
 * arguments the inner elements of `CovarArgs`, followed by either the inner
 * elements of `TrainTuple` or `PredTuple`. The return type of the `operator()`
 * method should be a type inheriting from `Eigen::EigenBase` with dynamic sized
 *  rows and columns.
 * @tparam Eta A type for the overdispersion parameter.
 * @tparam ThetaMatrix A type inheriting from `Eigen::EigenBase` with dynamic
 * sized rows and 1 column.
 * @tparam RNG A valid boost rng type
 * @tparam TrainTuple A tuple of types to passed as the end arguments of
 * `CovarFun::operator()`
 * @tparam PredTuple  A tuple of types to passed as the end arguments of
 * `CovarFun::operator()`
 * @tparam CovarArgs A tuple of types to passed as the first arguments of
 * `CovarFun::operator()`
 * @param y Observed counts.
 * @param y_index Index indicating which group each observation belongs to.
 * @param eta Overdisperison parameter.
 * @param theta_0 Initial guess for the Newton solver.
 * @param covariance_function Function that returns prior covariance matrix.
 * @param covar_args arguments for the covariance function.
 * @param train_tuple additional arguments for the covariance function,
 *                    e.g. covariates which correspond to observed data.
 * @param pred_tuple additional arguments for the covariance function,
 *                   e.g. covariates for out-of-sample data.
 * @param rng seed for rng.
 * @param msgs message stream for the covariance and likelihood function.
 */
template <typename CovarFun, typename Eta, typename ThetaMatrix, class RNG,
          typename TrainTuple, typename PredTuple, typename CovarArgs,
          require_eigen_t<ThetaMatrix>* = nullptr>
inline Eigen::VectorXd laplace_marginal_neg_binomial_2_log_rng(
    const std::vector<int>& y, const std::vector<int>& y_index, const Eta& eta,
    const ThetaMatrix& theta_0, CovarFun&& covariance_function,
    CovarArgs&& covar_args, TrainTuple&& train_tuple, PredTuple&& pred_tuple,
    RNG& rng, std::ostream* msgs) {
  return laplace_base_rng(neg_binomial_2_log_likelihood{},
                          std::forward_as_tuple(eta, y, y_index), theta_0,
                          std::forward<CovarFun>(covariance_function),
                          std::forward<CovarArgs>(covar_args),
                          std::forward<TrainTuple>(train_tuple),
                          std::forward<PredTuple>(pred_tuple), laplace_default_ops, rng, msgs);
}

}  // namespace math
}  // namespace stan

#endif
