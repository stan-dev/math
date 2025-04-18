#ifndef STAN_MATH_MIX_PROB_LAPLACE_LATENT_RNG_HPP
#define STAN_MATH_MIX_PROB_LAPLACE_LATENT_RNG_HPP

#include <stan/math/mix/functor/laplace_base_rng.hpp>
#include <stan/math/mix/functor/laplace_likelihood.hpp>
#include <stan/math/prim/fun/to_ref.hpp>

namespace stan {
namespace math {

/**
 * In a latent gaussian model,
 *
 *   theta ~ Normal(0, Sigma(phi))
 *   y ~ p(y|theta,phi)
 *
 * return a sample from the Laplace approximation to p(theta|y,phi),
 * where the log likelihood is given by L_f.
 * @tparam LLFunc Type of likelihood function.
 * @tparam LLArgs Type of arguments of likelihood function.
 * @tparam ThetaMatrix A type inheriting from `Eigen::EigenBase` with dynamic
 * sized rows and 1 column.
 * @tparam CovarFun A functor with an
 *  `operator()(CovarArgsElements..., {TrainTupleElements...|
 * PredTupleElements...})` method. The `operator()` method should accept as
 * arguments the inner elements of `CovarArgs`. The return type of the `operator()`
 * method should be a type inheriting from `Eigen::EigenBase` with dynamic sized
 *  rows and columns.
 * @tparam CovarArgs A tuple of types to passed as the first arguments of
 * `CovarFun::operator()`
 * @tparam RNG A valid boost rng type
 * @param L_f Function that returns log likelihood.
 * @param ll_args Arguments for likelihood function.
 * @param covariance_function Function that returns covariance function.
 * @param covar_args arguments for the covariance function.
 * @param theta_0 Initial guess for Newton solver.
 * @param tolerance Tolerated gradient norm for Newton solver.
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
template <typename LLFunc, typename LLArgs, typename ThetaMatrix,
          typename CovarFun, typename CovarArgs, typename RNG>
inline auto laplace_latent_tol_rng(
    LLFunc&& L_f, LLArgs&& ll_args, ThetaMatrix&& theta_0,
    CovarFun&& covariance_function, CovarArgs&& covar_args,
    const double tolerance,
    const int64_t max_num_steps, const int hessian_block_size, const int solver,
    const int max_steps_line_search, RNG& rng, std::ostream* msgs) {
  const laplace_options ops{hessian_block_size, solver, max_steps_line_search,
                            tolerance, max_num_steps};
  return laplace_base_rng(std::forward<LLFunc>(L_f),
                          std::forward<LLArgs>(ll_args),
                          std::forward<ThetaMatrix>(theta_0),
                          std::forward<CovarFun>(covariance_function),
                          std::forward<CovarArgs>(covar_args), ops, rng, msgs);
}

/**
 * In a latent gaussian model,
 *
 *   theta ~ Normal(0, Sigma(phi))
 *   y ~ p(y|theta,phi)
 *
 * return a sample from the Laplace approximation to p(theta|y,phi),
 * where the log likelihood is given by L_f.
 * @tparam LLFunc Type of likelihood function.
 * @tparam LLArgs Type of arguments of likelihood function.
 * @tparam ThetaMatrix A type inheriting from `Eigen::EigenBase` with dynamic
 * sized rows and 1 column.
 * @tparam CovarFun A functor with an
 *  `operator()(CovarArgsElements..., {TrainTupleElements...|
 * PredTupleElements...})` method. The `operator()` method should accept as
 * arguments the inner elements of `CovarArgs`, followed by either the inner
 * elements of `TrainTuple` or `PredTuple`. The return type of the `operator()`
 * method should be a type inheriting from `Eigen::EigenBase` with dynamic sized
 *  rows and columns.
 * @tparam RNG A valid boost rng type
 * @tparam CovarArgs A tuple of types to passed as the first arguments of
 * `CovarFun::operator()`
 * @param L_f Function that returns log likelihood.
 * @param ll_args Arguments for likelihood function.
 * @param covariance_function Function that returns covariance function.
 * @param covar_args arguments for the covariance function.
 * @param theta_0 Initial guess for Newton solver.
 * @param rng seed for rng.
 * @param msgs message stream for the covariance and likelihood function.
 */
template <typename LLFunc, typename LLArgs, typename ThetaMatrix,
          typename CovarFun, typename CovarArgs, typename RNG>
inline auto laplace_latent_rng(LLFunc&& L_f, LLArgs&& ll_args,
                                            ThetaMatrix&& theta_0,
                                            CovarFun&& covariance_function,
                                            CovarArgs&& covar_args, RNG& rng,
                                            std::ostream* msgs) {
  constexpr laplace_options ops{1, 1, 0, 1e-6, 100};
  return laplace_base_rng(std::forward<LLFunc>(L_f),
                          std::forward<LLArgs>(ll_args),
                          std::forward<ThetaMatrix>(theta_0),
                          std::forward<CovarFun>(covariance_function),
                          std::forward<CovarArgs>(covar_args), ops, rng, msgs);
}

}  // namespace math
}  // namespace stan

#endif
