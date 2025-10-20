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
 * @tparam ThetaVec A type inheriting from `Eigen::EigenBase`
 * with dynamic sized rows and 1 column.
 * @tparam Mean type of the mean of the latent normal distribution
 * \laplace_common_template_args
 * @tparam RNG A valid boost rng type
 * @param[in] y Observed counts.
 * @param[in] y_index Index indicating which group each observation belongs to.
 * @param[in] mean The mean of the latent normal variable.
 * \laplace_common_args
 * \laplace_options
 * \rng_arg
 * \msg_arg
 */
template <typename ThetaVec, typename Mean, typename CovarFun,
          typename CovarArgs, typename RNG,
          require_eigen_vector_t<ThetaVec>* = nullptr>
inline Eigen::VectorXd laplace_latent_tol_poisson_log_rng(
    const std::vector<int>& y, const std::vector<int>& y_index, Mean&& mean,
    CovarFun&& covariance_function, CovarArgs&& covar_args, ThetaVec&& theta_0,
    const double tolerance, const int max_num_steps,
    const int hessian_block_size, const int solver,
    const int max_steps_line_search, RNG& rng, std::ostream* msgs) {
  laplace_options_user_supplied ops{hessian_block_size,    solver,
                                    max_steps_line_search, tolerance,
                                    max_num_steps,         value_of(theta_0)};
  return laplace_base_rng(
      poisson_log_likelihood{},
      std::forward_as_tuple(y, y_index, std::forward<Mean>(mean)),
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
 * @tparam Mean type of the mean of the latent normal distribution
 * \laplace_common_template_args
 * @tparam RNG A valid boost rng type
 * @param[in] y Observed counts.
 * @param[in] y_index Index indicating which group each observation belongs to.
 * @param[in] mean The mean of the latent normal variable.
 * \laplace_common_args
 * \rng_arg
 * \msg_arg
 */
template <typename CovarFun, typename CovarArgs, typename RNG, typename Mean>
inline Eigen::VectorXd laplace_latent_poisson_log_rng(
    const std::vector<int>& y, const std::vector<int>& y_index, Mean&& mean,
    CovarFun&& covariance_function, CovarArgs&& covar_args, RNG& rng,
    std::ostream* msgs) {
  return laplace_base_rng(
      poisson_log_likelihood{},
      std::forward_as_tuple(y, y_index, std::forward<Mean>(mean)),
      std::forward<CovarFun>(covariance_function),
      std::forward<CovarArgs>(covar_args), laplace_options_default{}, rng,
      msgs);
}

}  // namespace math
}  // namespace stan

#endif
