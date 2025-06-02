#ifndef STAN_MATH_MIX_PROB_LAPLACE_LATENT_POISSON_LOG_RNG_HPP
#define STAN_MATH_MIX_PROB_LAPLACE_LATENT_POISSON_LOG_RNG_HPP

#include <stan/math/mix/functor/laplace_base_rng.hpp>
#include <stan/math/mix/functor/laplace_likelihood.hpp>
#include <stan/math/mix/prob/laplace_marginal_poisson_log_lpmf.hpp>

namespace stan {
namespace math {

/**
 * In a latent gaussian model,
 *
 *   theta ~ Normal(0, Sigma(phi))
 *   y ~ p(y|theta,phi)
 *
 * return a sample from the Laplace approximation to p(theta|y,phi).
 * The Laplace approximation is computed using a Newton solver.
 * In this specialized function, the likelihood p(y|theta) is a
 * \laplace_common_template_args
 * @tparam RNG A valid boost rng type
 * @param[in] y Observed counts.
 * @param[in] y_index Index indicating which group each observation belongs to.
 * \laplace_common_args
 * \laplace_options
 * \rng_arg
 * \msg_arg
 */
template <typename ThetaVec, typename CovarFun, typename CovarArgs,
          typename RNG, require_eigen_t<ThetaVec>* = nullptr>
inline Eigen::VectorXd laplace_latent_tol_poisson_log_rng(
    const std::vector<int>& y, const std::vector<int>& y_index,
    ThetaVec&& theta_0, CovarFun&& covariance_function, CovarArgs&& covar_args,
    const double tolerance, const int max_num_steps,
    const int hessian_block_size, const int solver,
    const int max_steps_line_search, RNG& rng, std::ostream* msgs) {
  laplace_options ops{hessian_block_size, solver, max_steps_line_search,
                      tolerance, max_num_steps};
  return laplace_base_rng(poisson_log_likelihood{},
                          std::forward_as_tuple(y, y_index),
                          std::forward<ThetaVec>(theta_0),
                          std::forward<CovarFun>(covariance_function),
                          std::forward<CovarArgs>(covar_args), ops, rng, msgs);
}

/**
 * In a latent gaussian model,
 *
 *   theta ~ Normal(0, Sigma(phi))
 *   y ~ p(y|theta,phi)
 *
 * return a sample from the Laplace approximation to p(theta|y,phi).
 * The Laplace approximation is computed using a Newton solver.
 * In this specialized function, the likelihood p(y|theta) is a
 * Poisson with a log link.
 * \laplace_common_template_args
 * @tparam RNG A valid boost rng type
 * @param[in] y Observed counts.
 * @param[in] y_index Index indicating which group each observation belongs to.
 * \laplace_common_args
 * \rng_arg
 * \msg_arg
 */
template <typename ThetaVec, typename CovarFun, typename CovarArgs,
          typename RNG, require_eigen_t<ThetaVec>* = nullptr>
inline Eigen::VectorXd laplace_latent_poisson_log_rng(
    const std::vector<int>& y, const std::vector<int>& y_index,
    const ThetaVec& theta_0, CovarFun&& covariance_function,
    CovarArgs&& covar_args, RNG& rng, std::ostream* msgs) {
  constexpr laplace_options ops{1, 1, 0, 1e-6, 100};
  return laplace_base_rng(poisson_log_likelihood{},
                          std::forward_as_tuple(y, y_index), theta_0,
                          std::forward<CovarFun>(covariance_function),
                          std::forward<CovarArgs>(covar_args), ops, rng, msgs);
}

}  // namespace math
}  // namespace stan

#endif
