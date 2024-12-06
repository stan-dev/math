#ifndef STAN_MATH_MIX_PROB_LAPLACE_RNG_HPP
#define STAN_MATH_MIX_PROB_LAPLACE_RNG_HPP

#include <stan/math/mix/functor/laplace_base_rng.hpp>
#include <stan/math/mix/functor/laplace_likelihood.hpp>
#include <stan/math/prim/fun/to_ref.hpp>

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
 * where the log likelihood is given by L_f.
 * @tparam LLFunc Type of likelihood function.
 * @tparam LLArgs Type of arguments of likelihood function.
 * @tparam ThetaMatrix A type inheriting from `Eigen::EigenBase` with dynamic
 * sized rows and 1 column.
 * @tparam CovarFun A functor with an
 *  `operator()(CovarArgsElements..., {TrainTupleElements...| PredTupleElements...})`
 *  method. The `operator()` method should accept as arguments the
 *  inner elements of `CovarArgs`, followed by either the inner elements of
 *  `TrainTuple` or `PredTuple`. The return type of the `operator()` method
 *  should be a type inheriting from `Eigen::EigenBase` with dynamic sized
 *  rows and columns.
 * @tparam RNG A valid boost rng type
 * @tparam CovarArgs A tuple of types to passed as the first arguments of `CovarFun::operator()`
 * @tparam TrainTuple A tuple of types to passed as the end arguments of `CovarFun::operator()`
 * @tparam PredTuple  A tuple of types to passed as the end arguments of `CovarFun::operator()`
 * @param L_f
 * @param l_args
 * @param covariance_function
 * @param theta_0
 * @param tolerance
 * @param max_num_steps
 * @param hessian_block_size
 * @param solver
 * @param max_steps_line_search
 * @param rng
 * @param msgs
 * @param train_tuple
 * @param pred_tuple
 * @param args
 */
template <typename LLFunc, typename LLArgs, typename ThetaMatrix, typename CovarFun,
          typename CovarArgs, typename TrainTuple, typename PredTuple, typename RNG>
inline Eigen::VectorXd laplace_marginal_tol_rng(
    LLFunc&& L_f, LLArgs&& l_args, const ThetaMatrix& theta_0, CovarFun&& covariance_function,
    CovarArgs&& covar_args, TrainTuple&& train_tuple, PredTuple&& pred_tuple,
    const double tolerance, const int64_t max_num_steps,
    const int hessian_block_size, const int solver,
    const int max_steps_line_search, RNG& rng, std::ostream* msgs) {
  const laplace_options ops{hessian_block_size, solver, max_steps_line_search,
                            tolerance, max_num_steps};
  return laplace_base_rng(std::forward<LLFunc>(L_f), std::forward<LLArgs>(l_args),
                          theta_0,
                          std::forward<CovarFun>(covariance_function),
                          std::forward<CovarArgs>(covar_args),
                          std::forward<TrainTuple>(train_tuple),
                          std::forward<PredTuple>(pred_tuple), ops, rng, msgs);
}

/**
 * In a latent gaussian model,
 *
 *   theta ~ Normal(theta | 0, Sigma(phi))
 *   y ~ pi(y | theta)
 *
 * return a multivariate normal random variate sampled
 * from the gaussian approximation of p(theta | y, phi)
 * where the log likelihood is given by L_f.
 * @tparam LLFunc
 * @tparam LLArgs
 * @tparam ThetaMatrix
 * @tparam CovarFun
 * @tparam RNG
 * @tparam TrainTuple
 * @tparam PredTuple
 * @tparam Args
 * @param L_f
 * @param l_args
 * @param covariance_function
 * @param theta_0
 * @param rng
 * @param msgs
 * @param train_tuple
 * @param pred_tuple
 * @param args
 */
template <typename LLFunc, typename LLArgs, typename ThetaMatrix, typename CovarFun,
          typename CovarArgs, typename TrainTuple, typename PredTuple, typename RNG>
inline Eigen::VectorXd laplace_marginal_rng(
    LLFunc&& L_f, LLArgs&& l_args, const ThetaMatrix& theta_0,
    CovarFun&& covariance_function, CovarArgs&& covar_args, TrainTuple&& train_tuple,
    PredTuple&& pred_tuple, RNG& rng, std::ostream* msgs) {
  constexpr laplace_options ops{1, 1, 0, 1e-6, 100};
  return laplace_base_rng(std::forward<LLFunc>(L_f), std::forward<LLArgs>(l_args),
                          theta_0,
                          std::forward<CovarFun>(covariance_function),
                          std::forward<CovarArgs>(covar_args),
                          std::forward<TrainTuple>(train_tuple),
                          std::forward<PredTuple>(pred_tuple),
                          ops, rng, msgs);
}


}  // namespace math
}  // namespace stan

#endif
