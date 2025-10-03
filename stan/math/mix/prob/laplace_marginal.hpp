#ifndef STAN_MATH_MIX_PROB_LAPLACE_MARGINAL_HPP
#define STAN_MATH_MIX_PROB_LAPLACE_MARGINAL_HPP

#include <stan/math/mix/functor/laplace_likelihood.hpp>
#include <stan/math/mix/functor/laplace_marginal_density.hpp>

namespace stan {
namespace math {
/**
 * Wrapper function around the laplace_marginal_density function.
 * Returns the marginal density p(y|phi) by marginalizing out
 * the latent gaussian variable theta, with a Laplace approximation.
 * See the laplace_marginal function for more details.
 * The data y is assumed to be real.
 * The function is "overloaded" below for the int y and lpmf case.
 *
 * @tparam propto If FALSE, log density is computed up to an additive const.
 * @tparam LFun The function which returns the log likelihood.
 * @tparam LArgs A tuple of arguments to the log likelihood.
 * @tparam EtaVec The type of the parameter arguments for the likelihood fn.
 * \laplace_common_template_args
 * @param[in] L_f a function which returns the log likelihood.
 * @param[in] l_args A tuple of arguments to pass to the log likelihood.
 * \laplace_common_args
 * \laplace_options
 * \msg_arg
 */
template <bool propto = false, typename LFun, typename LArgs, typename CovarFun,
          typename ThetaVec, typename CovarArgs,
          require_all_eigen_vector_t<ThetaVec>* = nullptr>
inline auto laplace_marginal_tol(
    LFun&& L_f, LArgs&& l_args, CovarFun&& covariance_function,
    CovarArgs&& covar_args, const ThetaVec& theta_0, double tolerance,
    int max_num_steps, const int hessian_block_size, const int solver,
    const int max_steps_line_search, std::ostream* msgs) {
  laplace_options_user_supplied ops{hessian_block_size,    solver,
                                    max_steps_line_search, tolerance,
                                    max_num_steps,         value_of(theta_0)};
  return laplace_marginal_density(
      std::forward<LFun>(L_f), std::forward<LArgs>(l_args),
      std::forward<CovarFun>(covariance_function),
      std::forward<CovarArgs>(covar_args), ops, msgs);
}

/**
 * Wrapper function around the laplace_marginal function.
 * Returns the marginal density p(y | phi) by marginalizing out
 * the latent gaussian variable, with a Laplace approximation.
 * See the laplace_marginal function for more details.
 * The data y is assumed to be real.
 * The function is "overloaded" below for the int y and lpmf case.
 *
 * @tparam propto If FALSE, log density is computed up to an additive const.
 * @tparam LFun The function which returns the log likelihood.
 * @tparam LArgs A tuple of arguments to the log likelihood.
 * \laplace_common_template_args
 * @param[in] L_f a function which returns the log likelihood.
 * @param[in] l_args A tuple of arguments to pass to the log likelihood
 * \laplace_common_args
 * \msg_arg
 */
template <bool propto = false, typename LFun, typename LArgs, typename CovarFun,
          typename CovarArgs>
inline auto laplace_marginal(LFun&& L_f, LArgs&& l_args,
                             CovarFun&& covariance_function,
                             CovarArgs&& covar_args, std::ostream* msgs) {
  return laplace_marginal_density(
      std::forward<LFun>(L_f), std::forward<LArgs>(l_args),
      std::forward<CovarFun>(covariance_function),
      std::forward<CovarArgs>(covar_args), laplace_options_default{}, msgs);
}

}  // namespace math
}  // namespace stan

#endif
