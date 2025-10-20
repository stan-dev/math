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
  template <typename ThetaVec, typename Eta, typename Mean,
            require_all_eigen_vector_t<ThetaVec>* = nullptr>
  inline auto operator()(const ThetaVec& theta, const Eta& eta,
                         const std::vector<int>& y,
                         const std::vector<int>& y_index, Mean&& mean,
                         std::ostream* pstream) const {
    Eigen::VectorXi n_per_group = Eigen::VectorXi::Zero(theta.size());
    Eigen::VectorXi counts_per_group = Eigen::VectorXi::Zero(theta.size());

    for (int i = 0; i < y.size(); i++) {
      n_per_group[y_index[i] - 1]++;
      counts_per_group[y_index[i] - 1] += y[i];
    }
    Eigen::Map<const Eigen::VectorXi> y_map(y.data(), y.size());

    auto theta_offset = to_ref(add(theta, mean));
    auto log_eta_plus_exp_theta = eval(log(add(eta, exp(theta_offset))));
    return sum(binomial_coefficient_log(subtract(add(y_map, eta), 1), y_map))
           + sum(
               add(elt_multiply(counts_per_group,
                                subtract(theta_offset, log_eta_plus_exp_theta)),
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
 * @tparam ThetaVec A type inheriting from `Eigen::EigenBase`
 * with dynamic sized rows and 1 column.
 * @tparam Mean type of the mean of the latent normal distribution
 * \laplace_common_template_args
 * @param[in] y observed counts.
 * @param[in] y_index group to which each observation belongs. Each group
 *            is parameterized by one element of theta.
 * @param[in] eta non-marginalized model parameters for the likelihood.
 * @param[in] mean the mean of the latent normal variable
 * \laplace_common_args
 * \laplace_options
 * \msg_arg
 */
template <bool propto = false, typename Eta, typename ThetaVec, typename Mean,
          typename CovarFun, typename CovarArgs,
          require_eigen_vector_t<ThetaVec>* = nullptr>
inline auto laplace_marginal_tol_neg_binomial_2_log_lpmf(
    const std::vector<int>& y, const std::vector<int>& y_index, const Eta& eta,
    Mean&& mean, CovarFun&& covariance_function, CovarArgs&& covar_args,
    const ThetaVec& theta_0, double tolerance, int max_num_steps,
    const int hessian_block_size, const int solver,
    const int max_steps_line_search, std::ostream* msgs) {
  laplace_options_user_supplied ops{hessian_block_size,    solver,
                                    max_steps_line_search, tolerance,
                                    max_num_steps,         value_of(theta_0)};
  return laplace_marginal_density(
      neg_binomial_2_log_likelihood{},
      std::forward_as_tuple(eta, y, y_index, std::forward<Mean>(mean)),
      std::forward<CovarFun>(covariance_function),
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
 * \laplace_common_template_args
 * @tparam Mean type of the mean of the latent normal distribution
 * @param[in] y observed counts.
 * @param[in] y_index group to which each observation belongs. Each group
 *            is parameterized by one element of theta.
 * @param[in] eta Parameter argument for likelihood function.
 * @param[in] mean the mean of the latent normal variable
 * \laplace_common_args
 * \msg_arg
 */
template <bool propto = false, typename Eta, typename Mean, typename CovarFun,
          typename CovarArgs>
inline auto laplace_marginal_neg_binomial_2_log_lpmf(
    const std::vector<int>& y, const std::vector<int>& y_index, const Eta& eta,
    Mean&& mean, CovarFun&& covariance_function, CovarArgs&& covar_args,
    std::ostream* msgs) {
  return laplace_marginal_density(
      neg_binomial_2_log_likelihood{},
      std::forward_as_tuple(eta, y, y_index, std::forward<Mean>(mean)),
      std::forward<CovarFun>(covariance_function),
      std::forward<CovarArgs>(covar_args), laplace_options_default{}, msgs);
}

struct neg_binomial_2_log_likelihood_summary {
  template <typename ThetaVec, typename Eta, typename Mean,
            require_eigen_vector_t<ThetaVec>* = nullptr>
  inline auto operator()(const ThetaVec& theta, const Eta& eta,
                         const std::vector<int>& y,
                         const std::vector<int>& n_per_group,
                         const std::vector<int>& counts_per_group, Mean&& mean,
                         std::ostream* pstream) const {
    Eigen::Map<const Eigen::VectorXi> y_map(y.data(), y.size());
    Eigen::Map<const Eigen::VectorXi> n_per_group_map(n_per_group.data(),
                                                      n_per_group.size());
    Eigen::Map<const Eigen::VectorXi> counts_per_group_map(
        counts_per_group.data(), counts_per_group.size());

    auto theta_offset = to_ref(add(theta, mean));
    auto log_eta_plus_exp_theta = eval(log(add(eta, exp(theta_offset))));

    return sum(binomial_coefficient_log(subtract(add(y_map, eta), 1.0), y_map))
           + sum(
               add(elt_multiply(counts_per_group_map,
                                subtract(theta_offset, log_eta_plus_exp_theta)),
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
 * @tparam ThetaVec A type inheriting from `Eigen::EigenBase`
 * with dynamic sized rows and 1 column.
 * @tparam Mean type of the mean of the latent normal distribution
 * \laplace_common_template_args
 * @param[in] y observations.
 * @param[in] n_per_group number of samples per group
 * @param[in] counts_per_group total counts per group
 * @param[in] eta non-marginalized model parameters for the likelihood.
 * @param[in] mean the mean of the latent normal variable
 * \laplace_common_args
 * \laplace_options
 * \msg_arg
 */
template <bool propto = false, typename Eta, typename ThetaVec, typename Mean,
          typename CovarFun, typename CovarArgs,
          require_eigen_vector_t<ThetaVec>* = nullptr>
inline auto laplace_marginal_tol_neg_binomial_2_log_summary_lpmf(
    const std::vector<int>& y, const std::vector<int>& n_per_group,
    const std::vector<int>& counts_per_group, const Eta& eta, Mean&& mean,
    CovarFun&& covariance_function, CovarArgs&& covar_args,
    const ThetaVec& theta_0, double tolerance, int max_num_steps,
    const int hessian_block_size, const int solver,
    const int max_steps_line_search, std::ostream* msgs) {
  laplace_options_user_supplied ops{hessian_block_size,    solver,
                                    max_steps_line_search, tolerance,
                                    max_num_steps,         value_of(theta_0)};
  return laplace_marginal_density(
      neg_binomial_2_log_likelihood_summary{},
      std::forward_as_tuple(eta, y, n_per_group, counts_per_group,
                            std::forward<Mean>(mean)),
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
 * @tparam Mean type of the mean of the latent normal distribution
 * \laplace_common_template_args
 * @param[in] y observations.
 * @param[in] n_per_group number of samples per group
 * @param[in] counts_per_group total counts per group
 * @param[in] eta non-marginalized model parameters for the likelihood.
 * @param[in] mean the mean of the latent normal variable
 * \laplace_common_args
 * \msg_arg
 */
template <bool propto = false, typename Eta, typename Mean, typename CovarFun,
          typename CovarArgs>
inline auto laplace_marginal_neg_binomial_2_log_summary_lpmf(
    const std::vector<int>& y, const std::vector<int>& n_per_group,
    const std::vector<int>& counts_per_group, const Eta& eta, Mean&& mean,
    CovarFun&& covariance_function, CovarArgs&& covar_args,
    std::ostream* msgs) {
  return laplace_marginal_density(
      neg_binomial_2_log_likelihood_summary{},
      std::forward_as_tuple(eta, y, n_per_group, counts_per_group,
                            std::forward<Mean>(mean)),
      std::forward<CovarFun>(covariance_function),
      std::forward<CovarArgs>(covar_args), laplace_options_default{}, msgs);
}

}  // namespace math
}  // namespace stan

#endif
