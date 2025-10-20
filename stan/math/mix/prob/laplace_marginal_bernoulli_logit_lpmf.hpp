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
  template <typename ThetaVec, typename YVec, typename Mean>
  inline auto operator()(const ThetaVec& theta, const YVec& y,
                         const std::vector<int>& delta_int, Mean&& mean,
                         std::ostream* pstream) const {
    auto theta_offset = to_ref(add(theta, mean));
    return sum(
        elt_multiply(theta_offset, y)
        - elt_multiply(to_vector(delta_int), log(add(1.0, exp(theta_offset)))));
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
 * @tparam ThetaVec A type inheriting from `Eigen::EigenBase`
 * with dynamic sized rows and 1 column.
 * @tparam Mean type of the mean of the latent normal distribution
 * \laplace_common_template_args
 * @param[in] y total counts per group. Second sufficient statistics.
 * @param[in] n_samples number of samples per group. First sufficient
 * statistics.
 * @param[in] mean the mean of the latent normal variable.
 * \laplace_common_args
 * \laplace_options
 * \msg_arg
 */
template <bool propto = false, typename ThetaVec, typename Mean,
          typename CovarFun, typename CovarArgs,
          require_eigen_vector_t<ThetaVec>* = nullptr>
inline auto laplace_marginal_tol_bernoulli_logit_lpmf(
    const std::vector<int>& y, const std::vector<int>& n_samples, Mean&& mean,
    CovarFun&& covariance_function, CovarArgs&& covar_args,
    const ThetaVec& theta_0, double tolerance, int max_num_steps,
    const int hessian_block_size, const int solver,
    const int max_steps_line_search, std::ostream* msgs) {
  laplace_options_user_supplied ops{hessian_block_size,    solver,
                                    max_steps_line_search, tolerance,
                                    max_num_steps,         value_of(theta_0)};
  return laplace_marginal_density(
      bernoulli_logit_likelihood{},
      std::forward_as_tuple(to_vector(y), n_samples, std::forward<Mean>(mean)),
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
 * @tparam Mean type of the mean of the latent normal distribution
 * \laplace_common_template_args
 * @param[in] y total counts per group. Second sufficient statistics.
 * @param[in] n_samples number of samples per group. First sufficient
 * statistics.
 * @param[in] mean the mean of the latent normal variable.
 * \laplace_common_args
 * \msg_arg
 */
template <bool propto = false, typename Mean, typename CovarFun,
          typename CovarArgs>
inline auto laplace_marginal_bernoulli_logit_lpmf(
    const std::vector<int>& y, const std::vector<int>& n_samples, Mean&& mean,
    CovarFun&& covariance_function, CovarArgs&& covar_args,
    std::ostream* msgs) {
  return laplace_marginal_density(
      bernoulli_logit_likelihood{},
      std::forward_as_tuple(to_vector(y), n_samples, std::forward<Mean>(mean)),
      std::forward<CovarFun>(covariance_function),
      std::forward<CovarArgs>(covar_args), laplace_options_default{}, msgs);
}

}  // namespace math
}  // namespace stan

#endif
