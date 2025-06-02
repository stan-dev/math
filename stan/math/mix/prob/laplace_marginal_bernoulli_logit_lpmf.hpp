#ifndef STAN_MATH_MIX_PROB_LAPLACE_MARGINAL_BERNOULLI_LOGIT_LPMF_HPP
#define STAN_MATH_MIX_PROB_LAPLACE_MARGINAL_BERNOULLI_LOGIT_LPMF_HPP

#include <stan/math/mix/functor/laplace_likelihood.hpp>
#include <stan/math/mix/functor/laplace_marginal_density.hpp>
#include <stan/math/rev/core/operator_addition.hpp>
#include <stan/math/rev/core/operator_subtraction.hpp>
#include <stan/math/rev/fun/elt_multiply.hpp>
#include <stan/math/rev/fun/exp.hpp>
#include <stan/math/rev/fun/log.hpp>
#include <stan/math/rev/fun/multiply.hpp>
#include <stan/math/rev/fun/sum.hpp>
#include <stan/math/fwd/fun/exp.hpp>
#include <stan/math/fwd/fun/lgamma.hpp>
#include <stan/math/fwd/fun/log.hpp>
#include <stan/math/fwd/fun/sum.hpp>
#include <stan/math/prim/fun/binomial_coefficient_log.hpp>

namespace stan {
namespace math {

struct bernoulli_logit_likelihood {
  template <typename T_theta, typename YVec>
  inline auto operator()(const T_theta& theta, const YVec& y,
                         const std::vector<int>& delta_int,
                         std::ostream* pstream) const {
    return sum(elt_multiply(theta, y)
               - elt_multiply(to_vector(delta_int), log(add(1.0, exp(theta)))));
  }
};

/**
 * Wrapper function around the laplace_marginal function for
 * a logistic Bernoulli likelihood. Returns the marginal density
 * p(y | phi) by marginalizing out the latent gaussian variable,
 * with a Laplace approximation. See the laplace_marginal function
 * for more details.
 *
 * @tparam propto boolean ignored
 * \laplace_common_template_args
 * @param[in] y total counts per group. Second sufficient statistics.
 * @param[in] n_samples number of samples per group. First sufficient
 *            statistics.
 * \laplace_common_args
 * \laplace_options
 * \msg_arg
 */
template <bool propto = false, typename ThetaVec, typename CovarFun,
          typename CovarArgs, require_eigen_t<ThetaVec>* = nullptr>
inline auto laplace_marginal_tol_bernoulli_logit_lpmf(
    const std::vector<int>& y, const std::vector<int>& n_samples,
    const ThetaVec& theta_0, CovarFun&& covariance_function,
    CovarArgs&& covar_args, double tolerance, int max_num_steps,
    const int hessian_block_size, const int solver,
    const int max_steps_line_search, std::ostream* msgs) {
  laplace_options ops{hessian_block_size, solver, max_steps_line_search,
                      tolerance, max_num_steps};
  return laplace_marginal_density(
      bernoulli_logit_likelihood{},
      std::forward_as_tuple(to_vector(y), n_samples), theta_0,
      std::forward<CovarFun>(covariance_function),
      std::forward<CovarArgs>(covar_args), ops, msgs);
}

/**
 * Wrapper function around the laplace_marginal function for
 * a logistic Bernoulli likelihood. Returns the marginal density
 * p(y | phi) by marginalizing out the latent gaussian variable,
 * with a Laplace approximation. See the laplace_marginal function
 * for more details.
 *
 * @tparam propto boolean ignored
 * \laplace_common_template_args
 * @param[in] y total counts per group. Second sufficient statistics.
 * @param[in] n_samples number of samples per group. First sufficient
 *            statistics.
 * \laplace_common_args
 * \msg_arg
 */
template <bool propto = false, typename CovarFun, typename ThetaVec,
          typename CovarArgs, require_eigen_t<ThetaVec>* = nullptr>
inline auto laplace_marginal_bernoulli_logit_lpmf(
    const std::vector<int>& y, const std::vector<int>& n_samples,
    const ThetaVec& theta_0, CovarFun&& covariance_function,
    CovarArgs&& covar_args, std::ostream* msgs) {
  constexpr laplace_options ops{1, 1, 0, 1e-6, 100};
  return laplace_marginal_density(
      bernoulli_logit_likelihood{},
      std::forward_as_tuple(to_vector(y), n_samples), theta_0,
      std::forward<CovarFun>(covariance_function),
      std::forward<CovarArgs>(covar_args), ops, msgs);
}

}  // namespace math
}  // namespace stan

#endif
