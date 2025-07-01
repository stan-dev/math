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
 * @tparam ThetaVec A type inheriting from `Eigen::EigenBase`
 * with dynamic sized rows and 1 column.
 * \laplace_common_template_args
 * @tparam RNG A valid boost rng type
 * @param[in] L_f Function that returns log likelihood.
 * @param[in] ll_args Arguments for likelihood function.
 * \laplace_common_args
 * \laplace_options
 * \rng_arg
 * \msg_arg
 */
template <typename LLFunc, typename LLArgs, typename ThetaVec,
          typename CovarFun, typename CovarArgs, typename RNG>
inline auto laplace_latent_tol_rng(
    LLFunc&& L_f, LLArgs&& ll_args, CovarFun&& covariance_function,
    CovarArgs&& covar_args, ThetaVec&& theta_0, const double tolerance,
    const int max_num_steps, const int hessian_block_size, const int solver,
    const int max_steps_line_search, RNG& rng, std::ostream* msgs) {
  const laplace_options_user_supplied ops{
      hessian_block_size, solver,        max_steps_line_search,
      tolerance,          max_num_steps, value_of(theta_0)};
  return laplace_base_rng(std::forward<LLFunc>(L_f),
                          std::forward<LLArgs>(ll_args),
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
 * \laplace_common_template_args
 * @tparam RNG A valid boost rng type
 * @param[in] L_f Function that returns log likelihood.
 * @param[in] ll_args Arguments for likelihood function.
 * \laplace_common_args
 * \rng_arg
 * \msg_arg
 */
template <typename LLFunc, typename LLArgs, typename CovarFun,
          typename CovarArgs, typename RNG>
inline auto laplace_latent_rng(LLFunc&& L_f, LLArgs&& ll_args,
                               CovarFun&& covariance_function,
                               CovarArgs&& covar_args, RNG& rng,
                               std::ostream* msgs) {
  return laplace_base_rng(std::forward<LLFunc>(L_f),
                          std::forward<LLArgs>(ll_args),
                          std::forward<CovarFun>(covariance_function),
                          std::forward<CovarArgs>(covar_args),
                          laplace_options_default{}, rng, msgs);
}

}  // namespace math
}  // namespace stan

#endif
