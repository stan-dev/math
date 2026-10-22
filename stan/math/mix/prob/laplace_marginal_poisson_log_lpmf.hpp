#ifndef STAN_MATH_MIX_PROB_LAPLACE_MARGINAL_POISSON_LOG_LPMF_HPP
#define STAN_MATH_MIX_PROB_LAPLACE_MARGINAL_POISSON_LOG_LPMF_HPP

#include <stan/math/mix/functor/laplace_likelihood.hpp>
#include <type_traits>
#include <stan/math/mix/functor/laplace_marginal_density.hpp>
#include <stan/math/mix/meta.hpp>
#include <stan/math/rev/fun/dot_product.hpp>
#include <stan/math/rev/fun/exp.hpp>
#include <stan/math/fwd/fun/exp.hpp>
#include <stan/math/rev/fun/sum.hpp>
#include <stan/math/fwd/fun/sum.hpp>
#include <stan/math/fwd/fun/lgamma.hpp>
#include <stan/math/rev/fun/lgamma.hpp>
#include <stan/math/prim/fun/add.hpp>

namespace stan {
namespace math {

struct poisson_log_likelihood {
  /**
   * Returns the lpmf for a Poisson with a log link across
   * multiple groups. No need to compute the log normalizing constant.
   * @tparam Theta A type inheriting from `Eigen::EigenBase` with dynamic
   * sized rows and 1 column.
   * @tparam YVec A vector type containing integers.
   * @tparam Mean type of the mean of the latent normal distribution
   * @param[in] theta log Poisson rate for each group.
   * @param[in] y observed counts
   * @param[in] y_index group to which each observation belongs
   * return lpmf for a Poisson with a log link.
   * @param[in] mean the mean of the latent normal variable
   * \msg_arg
   */
  template <typename Theta, typename YVec, typename Mean,
            require_eigen_vector_t<Theta>* = nullptr>
  inline auto operator()(const Theta& theta, const YVec& y,
                         const std::vector<int>& y_index, Mean&& mean,
                         std::ostream* /*pstream*/) const {
    Eigen::VectorXd counts_per_group = Eigen::VectorXd::Zero(theta.size());
    Eigen::VectorXd n_per_group = Eigen::VectorXd::Zero(theta.size());

    double norm_constant = 0;
    for (size_t i = 0; i < y_index.size(); i++) {
      counts_per_group(y_index[i] - 1) += y[i];
      n_per_group(y_index[i] - 1) += 1;
      norm_constant -= lgamma(y[i] + 1.0);
    }

    auto theta_offset = to_ref(add(theta, mean));

    return norm_constant + dot_product(theta_offset, counts_per_group)
           - dot_product(n_per_group, exp(theta_offset));
  }
};

namespace laplace_likelihood {
/**
 * Closed-form derivatives of `poisson_log_likelihood` with respect to
 * `theta` (see `analytic_derivatives`).  With `m_g = n_g exp(theta_g +
 * mean_g)`, where `n_g` is the number of observations of group `g`, the
 * gradient is `counts_g - m_g`, the negative Hessian is `m_g`, and the third
 * derivative is `-m_g`.
 */
template <>
struct analytic_derivatives<poisson_log_likelihood> : std::true_type {
  // sum of the counts and number of observations of each group
  template <typename YVec>
  static void group_sums(Eigen::Index theta_size, const YVec& y,
                         const std::vector<int>& y_index,
                         Eigen::VectorXd& counts, Eigen::VectorXd& n) {
    counts = Eigen::VectorXd::Zero(theta_size);
    n = Eigen::VectorXd::Zero(theta_size);
    for (size_t i = 0; i < y_index.size(); i++) {
      counts(y_index[i] - 1) += y[i];
      n(y_index[i] - 1) += 1;
    }
  }

  template <typename Theta, typename YVec, typename Mean>
  static void rates(const Theta& theta, const YVec& y,
                    const std::vector<int>& y_index, const Mean& mean,
                    Eigen::VectorXd& counts, Eigen::VectorXd& m) {
    Eigen::VectorXd n;
    group_sums(theta.size(), y, y_index, counts, n);
    Eigen::VectorXd theta_offset = add(theta, mean);
    m = n.cwiseProduct(theta_offset.array().exp().matrix());
  }

  template <typename Theta, typename YVec, typename Mean>
  static Eigen::VectorXd theta_grad(const Theta& theta, const YVec& y,
                                    const std::vector<int>& y_index,
                                    const Mean& mean) {
    Eigen::VectorXd counts, m;
    rates(theta, y, y_index, mean, counts, m);
    return counts - m;
  }

  template <typename Theta, typename YVec, typename Mean>
  static Eigen::VectorXd neg_hessian_diag(const Theta& theta, const YVec& y,
                                          const std::vector<int>& y_index,
                                          const Mean& mean) {
    Eigen::VectorXd counts, m;
    rates(theta, y, y_index, mean, counts, m);
    return m;
  }

  template <typename Theta, typename YVec, typename Mean>
  static Eigen::VectorXd third_diff(const Theta& theta, const YVec& y,
                                    const std::vector<int>& y_index,
                                    const Mean& mean) {
    Eigen::VectorXd counts, m;
    rates(theta, y, y_index, mean, counts, m);
    return -m;
  }
};
}  // namespace laplace_likelihood

/**
 * Wrapper function around the laplace_marginal function for
 * a log poisson likelihood. Returns the marginal density
 * p(y | phi) by marginalizing out the latent gaussian variable,
 * with a Laplace approximation. See the laplace_marginal function
 * for more details.
 *
 * @tparam propto ignored
 * @tparam ThetaVec A type inheriting from `Eigen::EigenBase`
 * with dynamic sized rows and 1 column.
 * @tparam Mean type of the mean of the latent normal distribution
 * \laplace_common_template_args
 * @param[in] y observed counts
 * @param[in] y_index group to which each observation belongs
 * @param[in] mean the mean of the latent normal variable
 * \laplace_common_args
 * @param[in] hessian_block_size Block size for the Hessian approximation with
 * respect to the latent gaussian variable theta.
 * \laplace_options
 * \msg_arg
 */
template <bool propto = false, typename Mean, typename CovarFun,
          typename CovarArgs, typename OpsTuple>
inline auto laplace_marginal_tol_poisson_log_lpmf(
    const std::vector<int>& y, const std::vector<int>& y_index, Mean&& mean,
    int hessian_block_size, CovarFun&& covariance_function,
    CovarArgs&& covar_args, OpsTuple&& ops, std::ostream* msgs) {
  auto options
      = internal::tuple_to_laplace_options(std::forward<OpsTuple>(ops));
  options.hessian_block_size = hessian_block_size;
  return laplace_marginal_density(
      poisson_log_likelihood{},
      std::forward_as_tuple(y, y_index, std::forward<Mean>(mean)),
      covariance_function, std::forward<CovarArgs>(covar_args),
      std::move(options), msgs);
}

/**
 * Wrapper function around the laplace_marginal function for
 * a log poisson likelihood. Returns the marginal density
 * p(y | phi) by marginalizing out the latent gaussian variable,
 * with a Laplace approximation. See the laplace_marginal function
 * for more details.
 *
 * @tparam propto ignored
 * @tparam Mean type of the mean of the latent normal distribution
 * \laplace_common_template_args
 * @param[in] y observed counts
 * @param[in] y_index group to which each observation belongs
 * @param[in] mean the mean of the latent normal variable
 * \laplace_common_args
 * @param[in] hessian_block_size Block size for the Hessian approximation with
 * respect to the latent gaussian variable theta.
 * \msg_arg
 */
template <bool propto = false, typename CovarFun, typename CovarArgs,
          typename Mean>
inline auto laplace_marginal_poisson_log_lpmf(
    const std::vector<int>& y, const std::vector<int>& y_index, Mean&& mean,
    int hessian_block_size, CovarFun&& covariance_function,
    CovarArgs&& covar_args, std::ostream* msgs) {
  auto options = laplace_options_default{hessian_block_size};
  return laplace_marginal_density(
      poisson_log_likelihood{},
      std::forward_as_tuple(y, y_index, std::forward<Mean>(mean)),
      covariance_function, std::forward<CovarArgs>(covar_args), options, msgs);
}

}  // namespace math
}  // namespace stan

#endif
