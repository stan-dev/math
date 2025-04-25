#ifndef STAN_MATH_MIX_PROB_LAPLACE_MARGINAL_NEG_BINOMIAL_2_LOG_LPMF_HPP
#define STAN_MATH_MIX_PROB_LAPLACE_MARGINAL_NEG_BINOMIAL_2_LOG_LPMF_HPP

#include <stan/math/mix/functor/laplace_likelihood.hpp>
#include <stan/math/mix/functor/laplace_marginal_density.hpp>

#include <stan/math/rev/core/operator_addition.hpp>
#include <stan/math/rev/core/operator_multiplication.hpp>
#include <stan/math/rev/core/operator_subtraction.hpp>
#include <stan/math/rev/fun/dot_product.hpp>
#include <stan/math/rev/fun/elt_multiply.hpp>
#include <stan/math/rev/fun/lgamma.hpp>
#include <stan/math/rev/fun/log.hpp>
#include <stan/math/rev/fun/exp.hpp>
#include <stan/math/rev/fun/multiply.hpp>
#include <stan/math/rev/fun/sum.hpp>
#include <stan/math/fwd/fun/exp.hpp>
#include <stan/math/fwd/fun/lgamma.hpp>
#include <stan/math/fwd/fun/log.hpp>
#include <stan/math/fwd/fun/sum.hpp>
#include <stan/math/prim/fun/binomial_coefficient_log.hpp>

namespace stan {
namespace math {

struct neg_binomial_2_log_likelihood {
  template <typename T_theta, typename T_eta>
  inline return_type_t<T_theta, T_eta> operator()(
      const T_theta& theta, const T_eta& eta,
      const std::vector<int>& y, const std::vector<int>& y_index,
      std::ostream* pstream) const {
    Eigen::VectorXi n_per_group = Eigen::VectorXi::Zero(theta.size());
    Eigen::VectorXi counts_per_group = Eigen::VectorXi::Zero(theta.size());

    for (int i = 0; i < y.size(); i++) {
      n_per_group[y_index[i]]++;
      counts_per_group[y_index[i]] += y[i];
    }
    Eigen::Map<const Eigen::VectorXi> y_map(y.data(), y.size());
    auto log_eta_plus_exp_theta = eval(log(add(eta, exp(theta))));
    return sum(binomial_coefficient_log(subtract(add(y_map, eta), 1), y_map))
           + sum(add(elt_multiply(counts_per_group,
                                  subtract(theta, log_eta_plus_exp_theta)),
                     elt_multiply(multiply(n_per_group, eta),
                                  subtract(log(eta), log_eta_plus_exp_theta))));
  }
};

/**
 * Wrapper function around the laplace_marginal function for
 * a negative binomial likelihood. Uses the 2nd parameterization.
 * Returns the marginal density p(y|phi) by marginalizing
 * out the latent gaussian variable, with a Laplace approximation.
 * See the laplace_marginal function for more details.
 *
 * @tparam Eta The type of parameter arguments for the likelihood function.
 * @tparam ThetaVec A type inheriting from `Eigen::EigenBase` with dynamic
 * sized rows and 1 column.
 * @tparam CovarFun A functor with an
 *  `operator()(CovarArgsElements...)` method. The `operator()` method should
 * accept as arguments the inner elements of `CovarArgs`. The return type of the
 * `operator()` method should be a type inheriting from `Eigen::EigenBase` with
 * dynamic sized rows and columns.
 * @tparam CovarArgs A tuple of types to passed as the first arguments of
 * `CovarFun::operator()`
 * @param[in] y observed counts.
 * @param[in] y_index group to which each observation belongs. Each group
 *            is parameterized by one element of theta.
 * @param[in] eta non-marginalized model parameters for the likelihood.
 * @param[in] theta_0 the initial guess for the Laplace approximation.
 * @param[in] covariance_function a function which returns the prior covariance.
 * @param[in] covar_args arguments for the covariance function.
 * @param[in] tolerance controls the convergence criterion when finding
 *            the mode in the Laplace approximation.
 * @param[in] max_num_steps maximum number of steps before the Newton solver
 *            breaks and returns an error.
 * @param[in] hessian_block_size the size of the block for a block-diagonal
 *              Hessian of the log likelihood. If 0, the Hessian is stored
 *              inside a vector. If the Hessian is dense, this should be the
 *              size of the Hessian.
 * @param[in] solver Type of Newton solver. Each corresponds to a distinct
 *               choice of B matrix (i.e. application SWM formula):
 *               1. computes square-root of negative Hessian.
 *               2. computes square-root of covariance matrix.
 *               3. computes no square-root and uses LU decomposition.
 * @param[in] max_steps_line_search Number of steps after which the algorithm
 *                          gives up on doing a linesearch. If 0, no linesearch.
 * @param[in,out] msgs message stream for the covariance and likelihood
 * function.
 */
template <bool propto = false, typename Eta, typename ThetaVec,
          typename CovarFun, typename CovarArgs,
          require_all_eigen_vector_t<ThetaVec>* = nullptr>
inline auto laplace_marginal_tol_neg_binomial_2_log_lpmf(
    const std::vector<int>& y, const std::vector<int>& y_index, const Eta& eta,
    const ThetaVec& theta_0, CovarFun&& covariance_function,
    CovarArgs&& covar_args, double tolerance, int64_t max_num_steps,
    const int hessian_block_size, const int solver,
    const int max_steps_line_search, std::ostream* msgs) {
  laplace_options ops{hessian_block_size, solver, max_steps_line_search,
                      tolerance, max_num_steps};
  return laplace_marginal_density(
      neg_binomial_2_log_likelihood{}, std::forward_as_tuple(eta, y, y_index),
      theta_0, std::forward<CovarFun>(covariance_function),
      std::forward<CovarArgs>(covar_args), ops, msgs);
}

/**
 * Wrapper function around the laplace_marginal function for
 * a negative binomial likelihood. Uses the 2nd parameterization.
 * Returns the marginal density p(y | phi) by marginalizing
 * out the latent gaussian variable, with a Laplace approximation.
 * See the laplace_marginal function for more details.
 *
 * @tparam Eta The type of parameter arguments for the likelihood function.
 * @tparam ThetaVec A type inheriting from `Eigen::EigenBase` with dynamic
 * sized rows and 1 column.
 * @tparam CovarFun A functor with an
 *  `operator()(CovarArgsElements...)` method. The `operator()` method should
 * accept as arguments the inner elements of `CovarArgs`. The return type of the
 * `operator()` method should be a type inheriting from `Eigen::EigenBase` with
 * dynamic sized rows and columns.
 * @tparam CovarArgs A tuple of types to passed as the first arguments of
 * `CovarFun::operator()`
 * @param[in] y observed counts.
 * @param[in] y_index group to which each observation belongs. Each group
 *            is parameterized by one element of theta.
 * @param[in] sums Total number of counts per group.
 * @param[in] eta Parameter argument for likelihood function.
 * @param[in] theta_0 the initial guess for the Laplace approximation.
 * @param[in] covariance_function a function which returns the prior covariance.
 * @param[in] covar_args arguments for the covariance function.
 * @param[in, out] msgs  message stream for the covariance and likelihood
 * function.
 */
template <bool propto = false, typename Eta, typename ThetaVec,
          typename CovarFun, typename CovarArgs>
inline auto laplace_marginal_neg_binomial_2_log_lpmf(
    const std::vector<int>& y, const std::vector<int>& y_index, const Eta& eta,
    const ThetaVec& theta_0, CovarFun&& covariance_function,
    CovarArgs&& covar_args, std::ostream* msgs) {
  constexpr laplace_options ops{1, 1, 0, 1e-6, 100};
  return laplace_marginal_density(
      neg_binomial_2_log_likelihood{}, std::forward_as_tuple(eta, y, y_index),
      theta_0, std::forward<CovarFun>(covariance_function),
      std::forward<CovarArgs>(covar_args), ops, msgs);
}

struct neg_binomial_2_log_likelihood_summary {
  template <typename T_theta, typename T_eta>
  inline return_type_t<T_theta, T_eta> operator()(
      const T_theta& theta, const T_eta& eta, const std::vector<int>& y,
      const std::vector<int>& n_per_group,
      const std::vector<int>& counts_per_group, std::ostream* pstream) const {
    Eigen::Map<const Eigen::VectorXi> y_map(y.data(), y.size());
    Eigen::Map<const Eigen::VectorXi> n_per_group_map(n_per_group.data(),
                                                      n_per_group.size());
    Eigen::Map<const Eigen::VectorXi> counts_per_group_map(
        counts_per_group.data(), counts_per_group.size());
    auto log_eta_plus_exp_theta = eval(log(add(eta, exp(theta))));
    return sum(binomial_coefficient_log(subtract(add(y_map, eta), 1.0), y_map))
           + sum(add(elt_multiply(counts_per_group_map,
                                  subtract(theta, log_eta_plus_exp_theta)),
                     elt_multiply(multiply(n_per_group_map, eta),
                                  subtract(log(eta), log_eta_plus_exp_theta))));
  }
};

/**
 * Wrapper function around the laplace_marginal function for
 * a negative binomial likelihood. Uses the 2nd parameterization.
 * Returns the marginal density p(y|phi) by marginalizing
 * out the latent gaussian variable, with a Laplace approximation.
 * See the laplace_marginal function for more details.
 *
 * @tparam Eta The type of parameter arguments for the likelihood function.
 * @tparam ThetaVec A type inheriting from `Eigen::EigenBase` with dynamic
 * sized rows and 1 column.
 * @tparam CovarFun A functor with an
 *  `operator()(CovarArgsElements...)` method. The `operator()` method should
 * accept as arguments the inner elements of `CovarArgs`. The return type of the
 * `operator()` method should be a type inheriting from `Eigen::EigenBase` with
 * dynamic sized rows and columns.
 * @tparam CovarArgs A tuple of types to passed as the first arguments of
 * `CovarFun::operator()`
 * @param[in] y observations.
 * @param[in] n_per_group number of samples per group
 * @param[in] counts_per_group total counts per group
 * @param[in] eta non-marginalized model parameters for the likelihood.
 * @param[in] theta_0 the initial guess for the Laplace approximation.
 * @param[in] covariance_function a function which returns the prior covariance.
 * @param[in] covar_args arguments for the covariance function.
 * @param[in] tolerance controls the convergence criterion when finding
 *            the mode in the Laplace approximation.
 * @param[in] max_num_steps maximum number of steps before the Newton solver
 *            breaks and returns an error.
 * @param[in] hessian_block_size the size of the block for a block-diagonal
 *              Hessian of the log likelihood. If 0, the Hessian is stored
 *              inside a vector. If the Hessian is dense, this should be the
 *              size of the Hessian.
 * @param[in] solver Type of Newton solver. Each corresponds to a distinct
 *               choice of B matrix (i.e. application SWM formula):
 *               1. computes square-root of negative Hessian.
 *               2. computes square-root of covariance matrix.
 *               3. computes no square-root and uses LU decomposition.
 * @param[in] max_steps_line_search Number of steps after which the algorithm
 *                          gives up on doing a linesearch. If 0, no linesearch.
 * @param[in, out] msgs message stream for the covariance and likelihood
 * function.
 */
template <bool propto = false, typename Eta, typename ThetaVec,
          typename CovarFun, typename CovarArgs,
          require_all_eigen_vector_t<ThetaVec>* = nullptr>
inline auto laplace_marginal_tol_neg_binomial_2_log_summary_lpmf(
    const std::vector<int>& y, const std::vector<int>& n_per_group,
    const std::vector<int>& counts_per_group, const Eta& eta,
    const ThetaVec& theta_0, CovarFun&& covariance_function,
    CovarArgs&& covar_args, double tolerance, int64_t max_num_steps,
    const int hessian_block_size, const int solver,
    const int max_steps_line_search, std::ostream* msgs) {
  laplace_options ops{hessian_block_size, solver, max_steps_line_search,
                      tolerance, max_num_steps};
  return laplace_marginal_density(
      neg_binomial_2_log_likelihood_summary{},
      std::forward_as_tuple(eta, y, n_per_group, counts_per_group), theta_0,
      std::forward<CovarFun>(covariance_function),
      std::forward<CovarArgs>(covar_args), ops, msgs);
}

/**
 * Wrapper function around the laplace_marginal function for
 * a negative binomial likelihood. Uses the 2nd parameterization.
 * Returns the marginal density p(y|phi) by marginalizing
 * out the latent gaussian variable, with a Laplace approximation.
 * See the laplace_marginal function for more details.
 *
 * @tparam Eta The type of parameter arguments for the likelihood function.
 * @tparam ThetaVec A type inheriting from `Eigen::EigenBase` with dynamic
 * sized rows and 1 column.
 * @tparam CovarFun A functor with an
 *  `operator()(CovarArgsElements...)` method. The `operator()` method should
 * accept as arguments the inner elements of `CovarArgs`. The return type of the
 * `operator()` method should be a type inheriting from `Eigen::EigenBase` with
 * dynamic sized rows and columns.
 * @tparam CovarArgs A tuple of types to passed as the first arguments of
 * `CovarFun::operator()`
 * @param[in] y observations.
 * @param[in] n_per_group number of samples per group
 * @param[in] counts_per_group total counts per group
 * @param[in] eta non-marginalized model parameters for the likelihood.
 * @param[in] theta_0 the initial guess for the Laplace approximation.
 * @param[in] covariance_function a function which returns the prior covariance.
 * @param[in] covar_args arguments for the covariance function.
 * @param[in, out] msgs message stream for the covariance and likelihood
 * function.
 */
template <bool propto = false, typename Eta, typename ThetaVec,
          typename CovarFun, typename CovarArgs,
          require_all_eigen_vector_t<ThetaVec>* = nullptr>
inline auto laplace_marginal_neg_binomial_2_log_summary_lpmf(
    const std::vector<int>& y, const std::vector<int>& n_per_group,
    const std::vector<int>& counts_per_group, const Eta& eta,
    const ThetaVec& theta_0, CovarFun&& covariance_function,
    CovarArgs&& covar_args, std::ostream* msgs) {
  constexpr laplace_options ops{1, 1, 0, 1e-6, 100};
  return laplace_marginal_density(
      neg_binomial_2_log_likelihood_summary{},
      std::forward_as_tuple(eta, y, n_per_group, counts_per_group), theta_0,
      std::forward<CovarFun>(covariance_function),
      std::forward<CovarArgs>(covar_args), ops, msgs);
}

}  // namespace math
}  // namespace stan

#endif
