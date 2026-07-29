#ifndef STAN_MATH_MIX_PROB_LAPLACE_LATENT_SOLVE_HPP
#define STAN_MATH_MIX_PROB_LAPLACE_LATENT_SOLVE_HPP

#include <stan/math/mix/functor/laplace_base_rng.hpp>

namespace stan {
namespace math {

namespace internal {
// Placeholder RNG type: used whenever laplace_base_rng's
// sampling branch is compiled out via `if constexpr`. In this case
// no rng is needed as no sampling is performed.
struct laplace_unused_rng {};

}  // namespace internal

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
 * \msg_arg
 */
template <typename LLFunc, typename LLArgs, typename CovarFun,
          typename CovarArgs, typename OpsTuple>
inline auto laplace_latent_solve_tol(LLFunc&& ll_fun, LLArgs&& ll_args,
                                     int hessian_block_size,
                                     CovarFun&& covariance_function,
                                     CovarArgs&& covar_args, OpsTuple&& ops,
                                     std::ostream* msgs) {
  auto options
      = internal::tuple_to_laplace_options(std::forward<OpsTuple>(ops));
  options.hessian_block_size = hessian_block_size;
  internal::laplace_unused_rng unused_rng;
  return laplace_base_rng<true>(std::forward<LLFunc>(ll_fun),
                                std::forward<LLArgs>(ll_args),
                                std::forward<CovarFun>(covariance_function),
                                std::forward<CovarArgs>(covar_args),
                                std::move(options), unused_rng, msgs);
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
 * \msg_arg
 */
template <typename LLFunc, typename LLArgs, typename CovarFun,
          typename CovarArgs>
inline auto laplace_latent_solve(LLFunc&& ll_fun, LLArgs&& ll_args,
                                 int hessian_block_size,
                                 CovarFun&& covariance_function,
                                 CovarArgs&& covar_args, std::ostream* msgs) {
  auto options = laplace_options_default{hessian_block_size};
  internal::laplace_unused_rng unused_rng;
  return laplace_base_rng<true>(std::forward<LLFunc>(ll_fun),
                                std::forward<LLArgs>(ll_args),
                                std::forward<CovarFun>(covariance_function),
                                std::forward<CovarArgs>(covar_args),
                                std::move(options), unused_rng, msgs);
}

}  // namespace math
}  // namespace stan

#endif
