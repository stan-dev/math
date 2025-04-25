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
 * @tparam ThetaVec The type of the initial guess, theta_0.
 * @tparam CovarFun A functor with an
 *  `operator()(CovarArgsElements...)` method. The `operator()` method should
 * accept as arguments the inner elements of `CovarArgs`. The return type of the
 * `operator()` method should be a type inheriting from `Eigen::EigenBase` with
 * dynamic sized rows and columns.
 * @tparam CovarArgs A tuple of types to passed as the first arguments of
 * `CovarFun::operator()`
 * @param[in] y total counts per group. Second sufficient statistics.
 * @param[in] n_samples number of samples per group. First sufficient
 *            statistics.
 * @param[in] theta_0 the initial guess for the Laplace approximation.
 * @param covariance_function Covariance function
 * @param covar_args arguments for the covariance function.
 * @param[in] tolerance controls the convergence criterion when finding
 *            the mode in the Laplace approximation.
 * @param[in] max_num_steps maximum number of steps before the Newton solver
 *            breaks and returns an error.
 * @param hessian_block_size Block size of Hessian of log likelihood w.r.t
 *                           latent Gaussian variable theta.
 * @param solver Type of Newton solver. Each corresponds to a distinct choice
 *                 of B matrix (i.e. application SWM formula):
 *               1. computes square-root of negative Hessian.
 *               2. computes square-root of covariance matrix.
 *               3. computes no square-root and uses LU decomposition.
 * @param max_steps_line_search Number of steps after which the algorithm gives
 *                              up on doing a linesearch. If 0, no linesearch.
 * @param msgs Rng number.
 * @param[in] args data for the covariance function.
 */
template <bool propto = false, typename CovarFun, typename ThetaVec,
          typename CovarArgs, require_eigen_t<ThetaVec>* = nullptr>
inline auto laplace_marginal_tol_bernoulli_logit_lpmf(
    const std::vector<int>& y, const std::vector<int>& n_samples,
    const ThetaVec& theta_0, CovarFun&& covariance_function,
    CovarArgs&& covar_args, double tolerance, int64_t max_num_steps,
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
 * @tparam CovarF Type of structure for covariance function.
 * @tparam ThetaVec The type of the initial guess, theta_0.
 * @tparam CovarFun A functor with an
 *  `operator()(CovarArgsElements...)` method. The `operator()` method should
 * accept as arguments the inner elements of `CovarArgs`. The return type of the
 * `operator()` method should be a type inheriting from `Eigen::EigenBase` with
 * dynamic sized rows and columns.
 * @tparam CovarArgs A tuple of types to passed as the first arguments of
 * `CovarFun::operator()`
 * @param[in] y total counts per group. Second sufficient statistics.
 * @param[in] n_samples number of samples per group. First sufficient
 *            statistics.
 * @param[in] theta_0 the initial guess for the Laplace approximation.
 * @param covariance_function
 * @param covar_args arguments for the covariance function.
 * @param msgs Streaming message for covariance functions.
 * @param[in] args data for the covariance function.
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
