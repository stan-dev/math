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
 * @tparam ThetaVec A type inheriting from `Eigen::EigenBase`
 * with dynamic sized rows and 1 column.
 * @tparam Mean type of the mean of the latent normal distribution
 * \laplace_common_template_args
 * @tparam RNG A valid boost rng type
 * @param[in] y Vector Vector of total number of trials with a positive outcome.
 * @param[in] n_samples Vector of number of trials.
 * @param[in] mean the mean of the latent normal variable.
 * \laplace_common_args
 * \laplace_options
 * \rng_arg
 * \msg_arg
 */
template <typename ThetaVec, typename Mean, typename CovarFun,
          typename CovarArgs, typename RNG,
          require_eigen_vector_t<ThetaVec>* = nullptr>
inline Eigen::VectorXd laplace_latent_tol_bernoulli_logit_rng(
    const std::vector<int>& y, const std::vector<int>& n_samples, Mean&& mean,
    CovarFun&& covariance_function, CovarArgs&& covar_args, ThetaVec&& theta_0,
    const double tolerance, const int max_num_steps,
    const int hessian_block_size, const int solver,
    const int max_steps_line_search, RNG& rng, std::ostream* msgs) {
  laplace_options_user_supplied ops{hessian_block_size,    solver,
                                    max_steps_line_search, tolerance,
                                    max_num_steps,         value_of(theta_0)};
  return laplace_base_rng(
      bernoulli_logit_likelihood{},
      std::forward_as_tuple(to_vector(y), n_samples, std::forward<Mean>(mean)),
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
 * @tparam Mean type of the mean of the latent normal distribution
 * \laplace_common_template_args
 * @tparam RNG A valid boost rng type
 * @param[in] y Vector Vector of total number of trials with a positive outcome.
 * @param[in] n_samples Vector of number of trials.
 * @param[in] mean the mean of the latent normal variable.
 * \laplace_common_args
 * \rng_arg
 * \msg_arg
 */
template <typename Mean, typename CovarFun, typename CovarArgs, typename RNG>
inline Eigen::VectorXd laplace_latent_bernoulli_logit_rng(
    const std::vector<int>& y, const std::vector<int>& n_samples, Mean&& mean,
    CovarFun&& covariance_function, CovarArgs&& covar_args, RNG& rng,
    std::ostream* msgs) {
  return laplace_base_rng(
      bernoulli_logit_likelihood{},
      std::forward_as_tuple(to_vector(y), n_samples, std::forward<Mean>(mean)),
      std::forward<CovarFun>(covariance_function),
      std::forward<CovarArgs>(covar_args), laplace_options_default{}, rng,
      msgs);
}

}  // namespace math
}  // namespace stan

#endif
