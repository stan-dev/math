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
 * @tparam Mean type of the mean of the latent normal distribution
 * \laplace_common_template_args
 * @tparam RNG A valid boost rng type
 * @param[in] y binary observations.
 * @param[in] y_index group to which each observation belongs.
 * @param[in] mean the mean of the latent normal variable.
 * \laplace_common_args
 * @param[in] hessian_block_size Block size for the Hessian approximation with
 * respect to the latent gaussian variable theta.
 * \laplace_options
 * \rng_arg
 * \msg_arg
 */
template <typename Mean, typename CovarFun, typename CovarArgs,
          typename OpsTuple, typename RNG>
inline Eigen::VectorXd laplace_latent_tol_bernoulli_logit_rng(
    const std::vector<int>& y, const std::vector<int>& y_index, Mean&& mean,
    int hessian_block_size, CovarFun&& covariance_function,
    CovarArgs&& covar_args, OpsTuple&& ops, RNG& rng, std::ostream* msgs) {
  auto options
      = internal::tuple_to_laplace_options(std::forward<OpsTuple>(ops));
  options.hessian_block_size = hessian_block_size;
  return laplace_base_rng(
      bernoulli_logit_likelihood{},
      std::forward_as_tuple(to_vector(y), y_index, std::forward<Mean>(mean)),
      std::forward<CovarFun>(covariance_function),
      std::forward<CovarArgs>(covar_args), std::move(options), rng, msgs);
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
 * @tparam Mean type of the mean of the latent normal distribution
 * \laplace_common_template_args
 * @tparam RNG A valid boost rng type
 * @param[in] y binary observations
 * @param[in] y_index group to which each observation belongs.
 * @param[in] mean the mean of the latent normal variable.
 * \laplace_common_args
 * @param[in] hessian_block_size Block size for the Hessian approximation with
 * respect to the latent gaussian variable theta.
 * \rng_arg
 * \msg_arg
 */
template <typename Mean, typename CovarFun, typename CovarArgs, typename RNG>
inline Eigen::VectorXd laplace_latent_bernoulli_logit_rng(
    const std::vector<int>& y, const std::vector<int>& y_index, Mean&& mean,
    int hessian_block_size, CovarFun&& covariance_function,
    CovarArgs&& covar_args, RNG& rng, std::ostream* msgs) {
  auto options = laplace_options_default{hessian_block_size};
  return laplace_base_rng(
      bernoulli_logit_likelihood{},
      std::forward_as_tuple(to_vector(y), y_index, std::forward<Mean>(mean)),
      std::forward<CovarFun>(covariance_function),
      std::forward<CovarArgs>(covar_args), options, rng, msgs);
}

}  // namespace math
}  // namespace stan

#endif
