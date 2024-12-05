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
 * @tparam LFun
 * @tparam LArgs
 * @tparam ThetaVec
 * @tparam CovarFun
 * @tparam RNG
 * @tparam TrainTuple
 * @tparam PredTuple
 * @tparam Args
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
template <typename LFun, typename LArgs, typename CovarFun, typename ThetaVec,
          typename RNG, typename TrainTuple, typename PredTuple,
          typename CovarArgs>
inline Eigen::VectorXd laplace_marginal_tol_rng(
    LFun&& L_f, LArgs&& l_args, const ThetaVec& theta_0, CovarFun&& covariance_function,
    CovarArgs&& covar_args, TrainTuple&& train_tuple, PredTuple&& pred_tuple,
    const double tolerance, const int64_t max_num_steps,
    const int hessian_block_size, const int solver,
    const int max_steps_line_search, RNG& rng, std::ostream* msgs) {
  const laplace_options ops{hessian_block_size, solver, max_steps_line_search,
                            tolerance, max_num_steps};
  return laplace_base_rng(std::forward<LFun>(L_f), std::forward<LArgs>(l_args),
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
 * @tparam LFun
 * @tparam LArgs
 * @tparam ThetaVec
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
template <typename LFun, typename LArgs, typename CovarFun,
          typename ThetaVec, typename RNG, typename TrainTuple,
          typename PredTuple, typename CovarArgs>
inline Eigen::VectorXd laplace_marginal_rng(
    LFun&& L_f, LArgs&& l_args, const ThetaVec& theta_0,
    CovarFun&& covariance_function, CovarArgs&& covar_args, TrainTuple&& train_tuple,
    PredTuple&& pred_tuple, RNG& rng, std::ostream* msgs) {
  constexpr laplace_options ops{1, 1, 0, 1e-6, 100};
  return laplace_base_rng(std::forward<LFun>(L_f), std::forward<LArgs>(l_args),
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
