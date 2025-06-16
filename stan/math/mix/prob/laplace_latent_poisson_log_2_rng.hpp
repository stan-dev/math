#ifndef STAN_MATH_MIX_PROB_LAPLACE_LATENT_POISSON_LOG_2_RNG_HPP
#define STAN_MATH_MIX_PROB_LAPLACE_LATENT_POISSON_LOG_2_RNG_HPP

#include <stan/math/mix/functor/laplace_base_rng.hpp>
#include <stan/math/mix/functor/laplace_likelihood.hpp>
#include <stan/math/mix/prob/laplace_marginal_poisson_log_2_lpmf.hpp>

namespace stan {
namespace math {

/**
 * Wrapper function around the laplace_marginal function for
 * a log poisson likelihood with exposure. Returns the marginal density
 * p(y | phi) by marginalizing out the latent gaussian variable,
 * with a Laplace approximation. See the laplace_marginal function
 * for more details.
 *
 * @tparam YeVec A type inheriting from `Eigen::EigenBase` with dynamic
 *  sized rows and 1 column.
 * \laplace_common_template_args
 * @tparam RNG A valid boost rng type
 * @param[in] y total counts per group. Second sufficient statistics.
 * @param[in] y_index group to which each observation belongs.
 * @param[in] ye the exposure for each group.
 * \laplace_common_args
 * \laplace_options
 * \rng_arg
 * \msg_arg
 */
template <typename YeVec, typename ThetaVec, typename CovarFun,
          typename CovarArgs, typename RNG,
          require_eigen_t<ThetaVec>* = nullptr>
inline auto laplace_latent_tol_poisson_2_log_rng(
    const std::vector<int>& y, const std::vector<int>& y_index, const YeVec& ye,
    CovarFun&& covariance_function, CovarArgs&& covar_args, ThetaVec&& theta_0,
    const double tolerance, const int max_num_steps,
    const int hessian_block_size, const int solver,
    const int max_steps_line_search, RNG& rng, std::ostream* msgs) {
  laplace_options ops{hessian_block_size, solver,        max_steps_line_search,
                      tolerance,          max_num_steps, value_of(theta_0)};
  return laplace_base_rng(poisson_log_2_likelihood{},
                          std::forward_as_tuple(y, y_index, ye),
                          std::forward<CovarFun>(covariance_function),
                          std::forward<CovarArgs>(covar_args), ops, rng, msgs);
}

/**
 * Wrapper function around the laplace_marginal function for
 * a log poisson likelihood with exposure. Returns the marginal density
 * p(y | phi) by marginalizing out the latent gaussian variable,
 * with a Laplace approximation. See the laplace_marginal function
 * for more details.
 *
 * @tparam YeVec A type inheriting from `Eigen::EigenBase` with dynamic
 *  sized rows and 1 column.
 * \laplace_common_template_args
 * @tparam RNG A valid boost rng type
 * @param[in] y total counts per group. Second sufficient statistics.
 * @param[in] y_index group to which each observation belongs.
 * @param[in] ye the exposure for each group.
 * \laplace_common_args
 * \rng_arg
 * \msg_arg
 */
template <typename YeVec, typename CovarFun, typename CovarArgs, typename RNG>
inline auto laplace_latent_poisson_2_log_rng(const std::vector<int>& y,
                                             const std::vector<int>& y_index,
                                             const YeVec& ye,
                                             CovarFun&& covariance_function,
                                             CovarArgs&& covar_args, RNG& rng,
                                             std::ostream* msgs) {
  constexpr laplace_options ops{1, 1, 0, 1e-6, 100};
  return laplace_base_rng(poisson_log_2_likelihood{},
                          std::forward_as_tuple(y, y_index, ye),
                          std::forward<CovarFun>(covariance_function),
                          std::forward<CovarArgs>(covar_args), ops, rng, msgs);
}

}  // namespace math
}  // namespace stan

#endif
