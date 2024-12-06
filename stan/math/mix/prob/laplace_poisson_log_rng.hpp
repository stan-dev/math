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
 * @tparam CovarFun A functor with an
 *  `operator()(CovarArgsElements..., {TrainTupleElements...| PredTupleElements...})`
 *  method. The `operator()` method should accept as arguments the
 *  inner elements of `CovarArgs`, followed by either the inner elements of
 *  `TrainTuple` or `PredTuple`. The return type of the `operator()` method
 *  should be a type inheriting from `Eigen::EigenBase` with dynamic sized
 *  rows and columns.
 * @tparam ThetaMatrix A type inheriting from `Eigen::EigenBase` with dynamic
 * sized rows and 1 column.
 * @tparam RNG A valid boost rng type
 * @tparam TrainTuple A tuple of types to passed as the end arguments of `CovarFun::operator()`
 * @tparam PredTuple  A tuple of types to passed as the end arguments of `CovarFun::operator()`
 * @tparam CovarArgs A tuple of types to passed as the first arguments of `CovarFun::operator()`
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
    CovarArgs&& covar_args, TrainTuple&& train_tuple, PredTuple&& pred_tuple,
    const double tolerance, const int64_t max_num_steps,
    const int hessian_block_size, const int solver,
    const int max_steps_line_search, RNG& rng, std::ostream* msgs) {
  laplace_options ops{hessian_block_size, solver, max_steps_line_search,
                      tolerance, max_num_steps};
  return laplace_base_rng(poisson_log_likelihood{},
                          std::forward_as_tuple(to_vector(y), n_samples),
                          theta_0,
                          std::forward<CovarFun>(covariance_function),
                          std::forward<CovarArgs>(covar_args),
                          std::forward<TrainTuple>(train_tuple),
                          std::forward<PredTuple>(pred_tuple),
                          ops, rng, msgs);
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
 * @tparam CovarFun A functor with an
 *  `operator()(CovarArgsElements..., {TrainTupleElements...| PredTupleElements...})`
 *  method. The `operator()` method should accept as arguments the
 *  inner elements of `CovarArgs`, followed by either the inner elements of
 *  `TrainTuple` or `PredTuple`. The return type of the `operator()` method
 *  should be a type inheriting from `Eigen::EigenBase` with dynamic sized
 *  rows and columns.
 * @tparam ThetaMatrix A type inheriting from `Eigen::EigenBase` with dynamic
 * sized rows and 1 column.
 * @tparam RNG A valid boost rng type
 * @tparam TrainTuple A tuple of types to passed as the end arguments of `CovarFun::operator()`
 * @tparam PredTuple  A tuple of types to passed as the end arguments of `CovarFun::operator()`
 * @tparam CovarArgs A tuple of types to passed as the first arguments of `CovarFun::operator()`
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
    CovarArgs&& covar_args, TrainTuple&& train_tuple, PredTuple&& pred_tuple,
    RNG& rng, std::ostream* msgs) {
  constexpr laplace_options ops{1, 1, 0, 1e-6, 100};
  return laplace_base_rng(poisson_log_likelihood{},
                          std::forward_as_tuple(to_vector(y), n_samples),
                          theta_0,
                          std::forward<CovarFun>(covariance_function),
                          std::forward<CovarArgs>(covar_args),
                          std::forward<TrainTuple>(train_tuple),
                          std::forward<PredTuple>(pred_tuple),
                          ops,rng, msgs);
}

}  // namespace math
}  // namespace stan

#endif
