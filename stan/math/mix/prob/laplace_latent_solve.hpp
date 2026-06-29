#ifndef STAN_MATH_MIX_PROB_LAPLACE_LATENT_SOLVE_HPP
#define STAN_MATH_MIX_PROB_LAPLACE_LATENT_SOLVE_HPP

#include <stan/math/mix/functor/laplace_base_rng.hpp>

namespace stan {
namespace math {

/**
 * In a latent gaussian model,
 *
 *   theta ~ Normal(0, Sigma(phi))
 *   y ~ p(y|theta,phi)
 *
 * returns the posterior mean and Cholesky factor from the Laplace
 * approximation to p(theta|y,phi), where the log likelihood is given by L_f.
 * @tparam LLFunc Type of likelihood function.
 * @tparam LLArgs Tuple of arguments types of likelihood function.
 * \laplace_common_template_args
 * @param ll_fun Likelihood function.
 * @param ll_args Arguments for likelihood function.
 * \laplace_common_args
 * @param[in] hessian_block_size Block size for the Hessian approximation with
 * respect to the latent gaussian variable theta.
 * \laplace_options
 * \rng_arg
 * \msg_arg
 */
template <typename LLFunc, typename LLArgs, typename CovarFun,
          typename CovarArgs, typename RNG, typename OpsTuple>
inline auto laplace_latent_tol_solve(LLFunc&& ll_fun, LLArgs&& ll_args,
                                     int hessian_block_size,
                                     CovarFun&& covariance_function,
                                     CovarArgs&& covar_args, OpsTuple&& ops,
                                     RNG& rng, std::ostream* msgs) {
  auto options
      = internal::tuple_to_laplace_options(std::forward<OpsTuple>(ops));
  options.hessian_block_size = hessian_block_size;
  return laplace_base_rng<true>(
      std::forward<LLFunc>(ll_fun), std::forward<LLArgs>(ll_args),
      std::forward<CovarFun>(covariance_function),
      std::forward<CovarArgs>(covar_args), std::move(options), rng, msgs);
}

/**
 * In a latent gaussian model,
 *
 *   theta ~ Normal(0, Sigma(phi))
 *   y ~ p(y|theta,phi)
 *
 * returns the posterior mean and Cholesky factor
 * from the Laplace approximation of p(theta | y, phi).
 * @tparam LLFunc Type of likelihood function.
 * @tparam LLArgs Tuple of arguments types of likelihood function.
 * \laplace_common_template_args
 * @tparam RNG A valid boost rng type
 * @param ll_fun Likelihood function.
 * @param ll_args Arguments for likelihood function.
 * \laplace_common_args
 * @param[in] hessian_block_size Block size for the Hessian approximation with
 * respect to the latent gaussian variable theta.
 * \rng_arg
 * \msg_arg
 */
template <typename LLFunc, typename LLArgs, typename CovarFun,
          typename CovarArgs, typename RNG>
inline auto laplace_latent_solve(LLFunc&& ll_fun, LLArgs&& ll_args,
                                 int hessian_block_size,
                                 CovarFun&& covariance_function,
                                 CovarArgs&& covar_args, RNG& rng,
                                 std::ostream* msgs) {
  auto options = laplace_options_default{hessian_block_size};
  return laplace_base_rng<true>(
      std::forward<LLFunc>(ll_fun), std::forward<LLArgs>(ll_args),
      std::forward<CovarFun>(covariance_function),
      std::forward<CovarArgs>(covar_args), std::move(options), rng, msgs);
}

}  // namespace math
}  // namespace stan

#endif
