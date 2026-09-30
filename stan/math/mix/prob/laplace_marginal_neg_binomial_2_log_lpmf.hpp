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
#include <stan/math/rev/fun/log_sum_exp.hpp>
#include <stan/math/rev/fun/exp.hpp>
#include <stan/math/rev/fun/multiply.hpp>
#include <stan/math/rev/fun/sum.hpp>
#include <stan/math/fwd/fun/exp.hpp>
#include <stan/math/fwd/fun/lgamma.hpp>
#include <stan/math/fwd/fun/log.hpp>
#include <stan/math/fwd/fun/log_sum_exp.hpp>
#include <stan/math/fwd/fun/sum.hpp>
#include <stan/math/prim/err/check_size_match.hpp>
#include <stan/math/prim/fun/binomial_coefficient_log.hpp>

namespace stan {
namespace math {

struct neg_binomial_2_log_likelihood {
  /**
   * Returns the lpmf for a negative binomial (2nd parameterization, log
   * link) across multiple groups. The dispersion `eta` is either a scalar
   * shared by all groups or a vector with one entry per group.
   */
  template <typename ThetaVec, typename Eta, typename Mean,
            require_all_eigen_vector_t<ThetaVec>* = nullptr>
  inline auto operator()(const ThetaVec& theta, const Eta& eta,
                         const std::vector<int>& y,
                         const std::vector<int>& y_index, Mean&& mean,
                         std::ostream* pstream) const {
    Eigen::VectorXi n_per_group = Eigen::VectorXi::Zero(theta.size());
    Eigen::VectorXi counts_per_group = Eigen::VectorXi::Zero(theta.size());

    for (size_t i = 0; i < y.size(); i++) {
      n_per_group[y_index[i] - 1]++;
      counts_per_group[y_index[i] - 1] += y[i];
    }
    Eigen::Map<const Eigen::VectorXi> y_map(y.data(), y.size());

    auto theta_offset = add(theta, mean);
    if constexpr (is_stan_scalar<Eta>::value) {
      // one dispersion shared by all groups
      auto log_eta = log(eta);
      auto lse = to_ref(log_sum_exp(theta_offset, log_eta));

      return sum(
                 binomial_coefficient_log(subtract(add(y_map, eta), 1.0),
                                          y_map))
             + sum(add(
                 // counts_per_group * (theta - log(eta + exp(theta)))
                 elt_multiply(counts_per_group, subtract(theta_offset, lse)),
                 // n_per_group * eta * (log(eta) - log(eta + exp(theta)))
                 elt_multiply(multiply(n_per_group, eta),
                              subtract(log_eta, lse))));
    } else {
      // one dispersion per group
      check_size_match("neg_binomial_2_log_likelihood", "eta", eta.size(),
                       "theta", theta.size());
      const auto& eta_ref = to_ref(eta);
      auto log_eta = to_ref(log(eta_ref));
      auto lse = to_ref(log_sum_exp(theta_offset, log_eta));
      // y + eta - 1 with the dispersion of each observation's group
      Eigen::Matrix<scalar_type_t<Eta>, Eigen::Dynamic, 1> y_plus_eta_m1(
          y.size());
      for (size_t i = 0; i < y.size(); ++i) {
        y_plus_eta_m1.coeffRef(i)
            = eta_ref.coeff(y_index[i] - 1) + (y[i] - 1.0);
      }

      return sum(binomial_coefficient_log(y_plus_eta_m1, y_map))
             + sum(add(
                 // counts_per_group * (theta - log(eta + exp(theta)))
                 elt_multiply(counts_per_group, subtract(theta_offset, lse)),
                 // n_per_group * eta * (log(eta) - log(eta + exp(theta)))
                 elt_multiply(elt_multiply(n_per_group, eta_ref),
                              subtract(log_eta, lse))));
    }
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
 * @param[in] eta the overdispersion parameter: a scalar shared by all
 *            groups, or a vector with one entry per group.
 * @param[in] mean the mean of the latent normal variable
 * \laplace_common_args
 * @param[in] hessian_block_size Block size for the Hessian approximation with
 * respect to the latent gaussian variable theta.
 * \laplace_options
 * \msg_arg
 */
template <bool propto = false, typename Eta, typename Mean, typename CovarFun,
          typename CovarArgs, typename OpsTuple>
inline auto laplace_marginal_tol_neg_binomial_2_log_lpmf(
    const std::vector<int>& y, const std::vector<int>& y_index, const Eta& eta,
    Mean&& mean, int hessian_block_size, CovarFun&& covariance_function,
    CovarArgs&& covar_args, OpsTuple&& ops, std::ostream* msgs) {
  auto options
      = internal::tuple_to_laplace_options(std::forward<OpsTuple>(ops));
  options.hessian_block_size = hessian_block_size;
  return laplace_marginal_density(
      neg_binomial_2_log_likelihood{},
      std::forward_as_tuple(eta, y, y_index, std::forward<Mean>(mean)),
      std::forward<CovarFun>(covariance_function),
      std::forward<CovarArgs>(covar_args), std::move(options), msgs);
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
 * @param[in] eta the overdispersion parameter: a scalar shared by all
 *            groups, or a vector with one entry per group.
 * @param[in] mean the mean of the latent normal variable
 * \laplace_common_args
 * @param[in] hessian_block_size Block size for the Hessian approximation with
 * respect to the latent gaussian variable theta.
 * \msg_arg
 */
template <bool propto = false, typename Eta, typename Mean, typename CovarFun,
          typename CovarArgs>
inline auto laplace_marginal_neg_binomial_2_log_lpmf(
    const std::vector<int>& y, const std::vector<int>& y_index, const Eta& eta,
    Mean&& mean, int hessian_block_size, CovarFun&& covariance_function,
    CovarArgs&& covar_args, std::ostream* msgs) {
  auto options = laplace_options_default{hessian_block_size};
  return laplace_marginal_density(
      neg_binomial_2_log_likelihood{},
      std::forward_as_tuple(eta, y, y_index, std::forward<Mean>(mean)),
      std::forward<CovarFun>(covariance_function),
      std::forward<CovarArgs>(covar_args), options, msgs);
}

struct neg_binomial_2_log_likelihood_summary {
  // Without a per-observation group index the dispersion cannot vary by
  // group, so `eta` must be a scalar here.
  template <typename ThetaVec, typename Eta, typename Mean,
            require_eigen_vector_t<ThetaVec>* = nullptr,
            require_stan_scalar_t<Eta>* = nullptr>
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

    auto theta_offset = add(theta, mean);
    auto log_eta = log(eta);
    auto lse = to_ref(log_sum_exp(theta_offset, log_eta));

    return sum(binomial_coefficient_log(subtract(add(y_map, eta), 1.0), y_map))
           + sum(add(
               // counts_per_group * (theta - log(eta + exp(theta)))
               elt_multiply(counts_per_group_map, subtract(theta_offset, lse)),
               // n_per_group * eta * (log(eta) - log(eta + exp(theta)))
               elt_multiply(multiply(n_per_group_map, eta),
                            subtract(log_eta, lse))));
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
 * @param[in] hessian_block_size Block size for the Hessian approximation with
 * respect to the latent gaussian variable theta.
 * \laplace_options
 * \msg_arg
 */
template <bool propto = false, typename Eta, typename Mean, typename CovarFun,
          typename CovarArgs, typename OpsTuple>
inline auto laplace_marginal_tol_neg_binomial_2_log_summary_lpmf(
    const std::vector<int>& y, const std::vector<int>& n_per_group,
    const std::vector<int>& counts_per_group, const Eta& eta, Mean&& mean,
    int hessian_block_size, CovarFun&& covariance_function,
    CovarArgs&& covar_args, OpsTuple&& ops, std::ostream* msgs) {
  auto options
      = internal::tuple_to_laplace_options(std::forward<OpsTuple>(ops));
  options.hessian_block_size = hessian_block_size;
  return laplace_marginal_density(
      neg_binomial_2_log_likelihood_summary{},
      std::forward_as_tuple(eta, y, n_per_group, counts_per_group,
                            std::forward<Mean>(mean)),
      std::forward<CovarFun>(covariance_function),
      std::forward<CovarArgs>(covar_args), std::move(options), msgs);
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
 * @param[in] hessian_block_size Block size for the Hessian approximation with
 * respect to the latent gaussian variable theta.
 * \msg_arg
 */
template <bool propto = false, typename Eta, typename Mean, typename CovarFun,
          typename CovarArgs>
inline auto laplace_marginal_neg_binomial_2_log_summary_lpmf(
    const std::vector<int>& y, const std::vector<int>& n_per_group,
    const std::vector<int>& counts_per_group, const Eta& eta, Mean&& mean,
    int hessian_block_size, CovarFun&& covariance_function,
    CovarArgs&& covar_args, std::ostream* msgs) {
  auto options = laplace_options_default{hessian_block_size};
  return laplace_marginal_density(
      neg_binomial_2_log_likelihood_summary{},
      std::forward_as_tuple(eta, y, n_per_group, counts_per_group,
                            std::forward<Mean>(mean)),
      std::forward<CovarFun>(covariance_function),
      std::forward<CovarArgs>(covar_args), options, msgs);
}

}  // namespace math
}  // namespace stan

#endif
