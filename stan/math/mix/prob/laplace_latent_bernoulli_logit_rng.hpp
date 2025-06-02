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
 * \laplace_common_template_args
 * @tparam RNG A valid boost rng type
 * @param[in] y Vector Vector of total number of trials with a positive outcome.
 * @param[in] n_samples Vector of number of trials.
 * \laplace_common_args
 * \laplace_options
 * \rng_arg
 * \msg_arg
 */
template <typename ThetaVec, typename CovarFun, typename CovarArgs,
          typename RNG, require_eigen_t<ThetaVec>* = nullptr>
inline Eigen::VectorXd laplace_latent_tol_bernoulli_logit_rng(
    const std::vector<int>& y, const std::vector<int>& n_samples,
    ThetaVec&& theta_0, CovarFun&& covariance_function, CovarArgs&& covar_args,
    const double tolerance, const int max_num_steps,
    const int hessian_block_size, const int solver,
    const int max_steps_line_search, RNG& rng, std::ostream* msgs) {
  laplace_options ops{hessian_block_size, solver, max_steps_line_search,
                      tolerance, max_num_steps};
  return laplace_base_rng(bernoulli_logit_likelihood{},
                          std::forward_as_tuple(to_vector(y), n_samples),
                          std::forward<ThetaVec>(theta_0),
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
 * \laplace_common_template_args
 * @tparam RNG A valid boost rng type
 * @param[in] y Vector Vector of total number of trials with a positive outcome.
 * @param[in] n_samples Vector of number of trials.
 * \laplace_common_args
 * \rng_arg
 * \msg_arg
 */
template <typename CovarFun, typename ThetaVec, typename CovarArgs,
          typename RNG, require_eigen_t<ThetaVec>* = nullptr>
inline Eigen::VectorXd laplace_latent_bernoulli_logit_rng(
    const std::vector<int>& y, const std::vector<int>& n_samples,
    ThetaVec&& theta_0, CovarFun&& covariance_function, CovarArgs&& covar_args,
    RNG& rng, std::ostream* msgs) {
  constexpr laplace_options ops{1, 1, 0, 1e-6, 100};
  return laplace_base_rng(bernoulli_logit_likelihood{},
                          std::forward_as_tuple(to_vector(y), n_samples),
                          std::forward<ThetaVec>(theta_0),
                          std::forward<CovarFun>(covariance_function),
                          std::forward<CovarArgs>(covar_args), ops, rng, msgs);
}

}  // namespace math
}  // namespace stan

#endif
