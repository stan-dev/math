#ifndef STAN_MATH_LAPLACE_LAPLACE_MARGINAL_NEG_BINOMIAL_2_HPP
#define STAN_MATH_LAPLACE_LAPLACE_MARGINAL_NEG_BINOMIAL_2_HPP

#include <stan/math/laplace/laplace_marginal.hpp>
//#include <stan/math/laplace/laplace_likelihood.hpp>
#include <stan/math/laplace/laplace_likelihood_general.hpp>

namespace stan {
namespace math {
/**
 * Wrapper function around the laplace_marginal function for
 * a negative binomial likelihood. Uses the 2nd parameterization.
 * Returns the marginal density p(y | phi) by marginalizing
 * out the latent gaussian variable, with a Laplace approximation.
 * See the laplace_marginal function for more details.
 *
 * @tparam T0 The type of the initial guess, theta_0.
 * @tparam T1 The type for the global parameter, phi.
 * @param[in] y observations.
 * @param[in] y_index group to which each observation belongs. Each group
 *            is parameterized by one element of theta.
 * @param[in] covariance a function which returns the prior covariance.
 * @param[in] phi model parameters for the covariance functor.
 * @param[in] eta non-marginalized model parameters for the likelihood.
 * @param[in] x data for the covariance functor.
 * @param[in] delta additional real data for the covariance functor.
 * @param[in] delta_int additional integer data for covariance functor.
 * @param[in] theta_0 the initial guess for the Laplace approximation.
 * @param[in] tolerance controls the convergence criterion when finding
 *            the mode in the Laplace approximation.
 * @param[in] max_num_steps maximum number of steps before the Newton solver
 *            breaks and returns an error.
 */
template <typename CovarFun, typename Eta, typename Theta0, typename... Args>
inline auto laplace_marginal_neg_binomial_2_log_lpmf(
    const std::vector<int>& y, const std::vector<int>& y_index,
    CovarFun&& covariance_function,
    const Eta& eta,
    const Theta0& theta_0,
    std::ostream* msgs = nullptr, double tolerance = 1e-6,
    long int max_num_steps = 100, const int hessian_block_size = 0,
    const int solver = 1, const int do_line_search = 0,
    const int max_steps_line_search = 10, Args&&... args) {
  return laplace_marginal_density(
      diff_neg_binomial_2_log(to_vector(y), y_index, theta_0.size()),
      std::forward<CovarFun>(covariance_function), eta,
      theta_0, msgs, tolerance, max_num_steps, hessian_block_size,
       solver, do_line_search, max_steps_line_search, std::forward<Args>(args)...);
}
}  // namespace math
}  // namespace stan

#endif
