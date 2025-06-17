#ifndef STAN_MATH_MIX_PROB_LAPLACE_MARGINAL_POISSON_LOG_LPMF_HPP
#define STAN_MATH_MIX_PROB_LAPLACE_MARGINAL_POISSON_LOG_LPMF_HPP

#include <stan/math/mix/functor/laplace_likelihood.hpp>
#include <stan/math/mix/functor/laplace_marginal_density.hpp>
#include <stan/math/rev/fun/dot_product.hpp>
#include <stan/math/rev/fun/exp.hpp>
#include <stan/math/rev/fun/lgamma.hpp>
#include <stan/math/rev/fun/sum.hpp>
#include <stan/math/fwd/fun/exp.hpp>
#include <stan/math/fwd/fun/lgamma.hpp>

namespace stan {
namespace math {

struct poisson_log_likelihood {
  /**
   * Returns the lpmf for a Poisson with a log link across
   * multiple groups. No need to compute the log normalizing constant.
   * @tparam T_theta Type of the log Poisson rate.
   * @tparam T_eta Type of the auxiliary parameter (not used here).
   * @param[in] theta log Poisson rate for each group.
   * @param[in] y observed counts
   * @param[in] y_index group to which each observation belongs
   * return lpmf for a Poisson with a log link.
   * @param[in] pstream
   */
  template <typename Theta, typename YVec,
            require_eigen_vector_t<Theta>* = nullptr>
  inline auto operator()(const Theta& theta, const YVec& y,
                         const std::vector<int>& y_index,
                         std::ostream* pstream) const {
    Eigen::VectorXd counts_per_group = Eigen::VectorXd::Zero(theta.size());
    Eigen::VectorXd n_per_group = Eigen::VectorXd::Zero(theta.size());

    for (int i = 0; i < theta.size(); i++) {
      counts_per_group(y_index[i]) += y[i];
      n_per_group(y_index[i]) += 1;
    }

    return -sum(lgamma(add(counts_per_group, 1)))
           + dot_product(theta, counts_per_group)
           - dot_product(n_per_group, exp(theta));
  }
};

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
 * \laplace_common_template_args
 * @param[in] y observed counts
 * @param[in] y_index group to which each observation belongs
 * \laplace_common_args
 * \laplace_options
 * \msg_arg
 */
template <bool propto = false, typename ThetaVec, typename CovarFun,
          typename CovarArgs, require_all_eigen_vector_t<ThetaVec>* = nullptr>
inline auto laplace_marginal_tol_poisson_log_lpmf(
    const std::vector<int>& y, const std::vector<int>& y_index,
    CovarFun&& covariance_function, CovarArgs&& covar_args,
    const ThetaVec& theta_0, double tolerance, int max_num_steps,
    const int hessian_block_size, const int solver,
    const int max_steps_line_search, std::ostream* msgs) {
  laplace_options ops{hessian_block_size, solver,        max_steps_line_search,
                      tolerance,          max_num_steps, value_of(theta_0)};
  return laplace_marginal_density(
      poisson_log_likelihood{}, std::forward_as_tuple(y, y_index),
      covariance_function, std::forward<CovarArgs>(covar_args), ops, msgs);
}

/**
 * Wrapper function around the laplace_marginal function for
 * a log poisson likelihood. Returns the marginal density
 * p(y | phi) by marginalizing out the latent gaussian variable,
 * with a Laplace approximation. See the laplace_marginal function
 * for more details.
 *
 * @tparam propto ignored
 * \laplace_common_template_args
 * @param[in] y observed counts
 * @param[in] y_index group to which each observation belongs
 * \laplace_common_args
 * \msg_arg
 */
template <bool propto = false, typename CovarFun, typename CovarArgs>
inline auto laplace_marginal_poisson_log_lpmf(const std::vector<int>& y,
                                              const std::vector<int>& y_index,
                                              CovarFun&& covariance_function,
                                              CovarArgs&& covar_args,
                                              std::ostream* msgs) {
  const laplace_options ops{1, 1, 0, 1e-6, 100, std::nullopt};
  return laplace_marginal_density(
      poisson_log_likelihood{}, std::forward_as_tuple(y, y_index),
      covariance_function, std::forward<CovarArgs>(covar_args), ops, msgs);
}

}  // namespace math
}  // namespace stan

#endif
