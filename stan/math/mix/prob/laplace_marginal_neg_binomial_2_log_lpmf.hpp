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

    auto theta_offset = add(theta, mean);
    auto log_eta = log(eta);
    auto lse = to_ref(log_sum_exp(theta_offset, log_eta));

    return sum(binomial_coefficient_log(subtract(add(y_map, eta), 1.0), y_map))
           + sum(add(
               // counts_per_group * (theta - log(eta + exp(theta)))
               elt_multiply(counts_per_group, subtract(theta_offset, lse)),
               // n_per_group * eta * (log(eta) - log(eta + exp(theta)))
               elt_multiply(multiply(n_per_group, eta),
                            subtract(log_eta, lse))));
  }

  /**
   * Compute gradient and negative Hessian of the neg_binomial_2_log
   * likelihood analytically, avoiding nested autodiff.
   *
   * @param theta Latent Gaussian variable (double).
   * @param hessian_block_size Size of each diagonal block (typically 1).
   * @param eta Dispersion parameter (scalar or 1-element vector).
   * @param y Observed counts.
   * @param y_index Group index for each observation.
   * @param mean Mean offset for theta.
   * @return pair of (gradient, negative Hessian) as (VectorXd, SparseMatrix).
   */
  template <typename Mean>
  inline auto diff(const Eigen::VectorXd& theta, int hessian_block_size,
                   const Eigen::VectorXd& eta, const std::vector<int>& y,
                   const std::vector<int>& y_index, Mean&& mean) const {
    const int theta_size = theta.size();
    const double eta_scalar = eta(0);

    Eigen::VectorXd sums = Eigen::VectorXd::Zero(theta_size);
    Eigen::VectorXd n_samples = Eigen::VectorXd::Zero(theta_size);
    for (size_t i = 0; i < y.size(); i++) {
      n_samples(y_index[i] - 1) += 1;
      sums(y_index[i] - 1) += y[i];
    }

    // theta + mean
    Eigen::VectorXd theta_offset = add(theta, value_of(mean));

    // exp(-theta_offset)
    Eigen::VectorXd exp_neg_theta = exp(-theta_offset);
    // sums + eta * n_samples
    Eigen::VectorXd sums_plus_n_eta = sums + eta_scalar * n_samples;
    // 1 + eta * exp(-theta_offset)
    Eigen::VectorXd one_plus_exp
        = Eigen::VectorXd::Ones(theta_size) + eta_scalar * exp_neg_theta;

    // Gradient: sums - (sums + eta * n) / (1 + eta * exp(-theta))
    Eigen::VectorXd gradient = sums - sums_plus_n_eta.cwiseQuotient(one_plus_exp);

    // Negative Hessian diagonal:
    //   eta * (sums + eta * n) * exp(-theta) / (1 + eta * exp(-theta))^2
    Eigen::VectorXd hessian_diag
        = eta_scalar
          * sums_plus_n_eta.cwiseProduct(
                exp_neg_theta.cwiseQuotient(one_plus_exp.cwiseProduct(one_plus_exp)));

    Eigen::SparseMatrix<double> hessian(theta_size, theta_size);
    hessian.reserve(
        Eigen::VectorXi::Constant(theta_size, hessian_block_size));
    for (int i = 0; i < theta_size; i++) {
      hessian.insert(i, i) = hessian_diag(i);
    }

    return std::make_pair(std::move(gradient), std::move(hessian));
  }

  /**
   * Compute the third derivative of the neg_binomial_2_log likelihood
   * w.r.t. theta analytically, avoiding fvar<fvar<var>> autodiff.
   *
   * The third derivative is:
   *   d^3/dtheta^3 log p(y|theta,eta) =
   *     -(sums + eta*n) * eta * exp(theta) * (eta - exp(theta))
   *     / (eta + exp(theta))^3
   *
   * @param theta Latent Gaussian variable (double).
   * @param eta Dispersion parameter (scalar or 1-element vector).
   * @param y Observed counts.
   * @param y_index Group index for each observation.
   * @param mean Mean offset for theta.
   * @return Third derivative as a VectorXd.
   */
  template <typename Mean>
  inline Eigen::VectorXd third_diff(const Eigen::VectorXd& theta,
                                    const Eigen::VectorXd& eta,
                                    const std::vector<int>& y,
                                    const std::vector<int>& y_index,
                                    Mean&& mean) const {
    const int theta_size = theta.size();
    const double eta_scalar = eta(0);

    Eigen::VectorXd sums = Eigen::VectorXd::Zero(theta_size);
    Eigen::VectorXd n_samples = Eigen::VectorXd::Zero(theta_size);
    for (size_t i = 0; i < y.size(); i++) {
      n_samples(y_index[i] - 1) += 1;
      sums(y_index[i] - 1) += y[i];
    }

    // theta + mean
    Eigen::VectorXd theta_offset = add(theta, value_of(mean));

    Eigen::VectorXd exp_theta = exp(theta_offset);
    Eigen::VectorXd eta_vec
        = Eigen::VectorXd::Constant(theta_size, eta_scalar);
    Eigen::VectorXd eta_plus_exp_theta = eta_vec + exp_theta;

    // -(sums + eta*n) * eta * exp(theta) * (eta - exp(theta))
    //   / (eta + exp(theta))^3
    Eigen::VectorXd eta_plus_exp_theta_sq
        = eta_plus_exp_theta.cwiseProduct(eta_plus_exp_theta);
    Eigen::VectorXd eta_plus_exp_theta_cubed
        = eta_plus_exp_theta_sq.cwiseProduct(eta_plus_exp_theta);

    return -((sums + eta_scalar * n_samples) * eta_scalar)
                .cwiseProduct(exp_theta.cwiseProduct(
                    (eta_vec - exp_theta)
                        .cwiseQuotient(eta_plus_exp_theta_cubed)));
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
 * @param[in] eta Parameter argument for likelihood function.
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
