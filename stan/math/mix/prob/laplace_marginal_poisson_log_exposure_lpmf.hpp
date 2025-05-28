#ifndef STAN_MATH_MIX_PROB_LAPLACE_MARGINAL_POISSON_LOG_EXPOSURE_LPMF_HPP
#define STAN_MATH_MIX_PROB_LAPLACE_MARGINAL_POISSON_LOG_EXPOSURE_LPMF_HPP

#include <stan/math/mix/functor/laplace_marginal_density.hpp>
#include <stan/math/mix/functor/laplace_likelihood.hpp>
#include <stan/math/rev/core/operator_addition.hpp>
#include <stan/math/rev/core/operator_multiplication.hpp>
#include <stan/math/rev/core/operator_subtraction.hpp>
#include <stan/math/rev/fun/dot_product.hpp>
#include <stan/math/rev/fun/exp.hpp>
#include <stan/math/rev/fun/lgamma.hpp>
#include <stan/math/rev/fun/log.hpp>
#include <stan/math/rev/fun/sum.hpp>
#include <stan/math/fwd/fun/exp.hpp>
#include <stan/math/fwd/fun/lgamma.hpp>
#include <stan/math/fwd/fun/sum.hpp>

namespace stan {
namespace math {

struct poisson_log_exposure_likelihood {
  /**
   * Returns the lpmf for a Poisson with a log link across
   * multiple groups. No need to compute the log normalizing constant.
   * Same as above, but includes a exposure term to correct the
   * log rate for each group.
   * @tparam Theta Type of the log Poisson rate.
   * @tparam Eta Type of the auxiliary parameter (not used here).
   * @param[in] theta log Poisson rate for each group.
   * @param[in] y First n elements contain the sum of counts in each group
   * @param[in] y_index group to which each observation belongs.
   * @param[in] ye next n elements the exposure in each group, where n is the
   * number of groups.
   * @param[in, out] pstream msgs that are not used here
   */
  template <typename Theta, typename YVec, typename YIndexVec, typename YeVec>
  inline auto operator()(const Theta& theta, const YVec& y,
                         const YIndexVec& y_index, const YeVec& ye,
                         std::ostream* /*pstream*/) const {
    Eigen::VectorXd y_vec = to_vector(y);
    Eigen::VectorXd counts_per_group = Eigen::VectorXd::Zero(theta.size());
    Eigen::VectorXd n_per_group = Eigen::VectorXd::Zero(theta.size());
    for (int i = 0; i < theta.size(); i++) {
      counts_per_group(y_index[i]) += y[i];
      n_per_group(y_index[i]) += 1;
    }
    // auto n_samples = to_vector(delta_int);
    auto shifted_mean = to_ref(add(theta, log(ye)));
    return -sum(lgamma(add(y_vec, 1))) + dot_product(shifted_mean, y_vec)
           - dot_product(n_per_group, exp(shifted_mean));
  }
};

/**
 * Wrapper function around the laplace_marginal function for
 * a log poisson likelihood. Returns the marginal density
 * p(y | phi) by marginalizing out the latent gaussian variable,
 * with a Laplace approximation. See the laplace_marginal function
 * for more details.
 *
 * @tparam propto boolean ignored
 * @tparam YeVec A type inheriting from `Eigen::EigenBase` with dynamic
 * sized rows and 1 column.
 * \laplace_common_template_args
 * @param[in] y total counts per group. Second sufficient statistics.
 * @param[in] y_index group to which each observation belongs.
 * @param[in] ye the exposure for each group.
 * \laplace_common_args
 * \laplace_options
 * \msg_arg
 */
template <bool propto = false, typename YeVec, typename ThetaVec,
          typename CovarFun, typename CovarArgs,
          require_all_eigen_vector_t<YeVec, ThetaVec>* = nullptr>
inline auto laplace_marginal_tol_poisson_2_log_lpmf(
    const std::vector<int>& y, const std::vector<int>& y_index, const YeVec& ye,
    const ThetaVec& theta_0, CovarFun&& covariance_function,
    CovarArgs&& covar_args, double tolerance, int64_t max_num_steps,
    const int hessian_block_size, const int solver,
    const int max_steps_line_search, std::ostream* msgs) {
  laplace_options ops{hessian_block_size, solver, max_steps_line_search,
                      tolerance, max_num_steps};
  return laplace_marginal_density(
      poisson_log_exposure_likelihood{}, std::forward_as_tuple(y, y_index, ye),
      theta_0, std::forward<CovarFun>(covariance_function),
      std::forward<CovarArgs>(covar_args), ops, msgs);
}

/**
 * Wrapper function around the laplace_marginal function for
 * a log poisson likelihood. Returns the marginal density
 * p(y | phi) by marginalizing out the latent gaussian variable,
 * with a Laplace approximation. See the laplace_marginal function
 * for more details.
 *
 * @tparam propto boolean ignored
 * @tparam YeVec The type for the global parameter, phi.
 * \laplace_common_template_args
 * @param[in] y total counts per group. Second sufficient statistics.
 * @param[in] y_index group to which each observation belongs.
 * @param[in] ye the exposure for each group.
 * \laplace_common_args
 * \msg_arg
 */
template <bool propto = false, typename YeVec, typename ThetaVec,
          typename CovarFun, typename CovarArgs,
          require_all_eigen_vector_t<YeVec, ThetaVec>* = nullptr>
inline auto laplace_marginal_poisson_2_log_lpmf(
    const std::vector<int>& y, const std::vector<int>& y_index, const YeVec& ye,
    const ThetaVec& theta_0, CovarFun&& covariance_function,
    CovarArgs&& covar_args, std::ostream* msgs) {
  constexpr laplace_options ops{1, 1, 0, 1e-6, 100};
  return laplace_marginal_density(
      poisson_log_exposure_likelihood{}, std::forward_as_tuple(y, y_index, ye),
      theta_0, std::forward<CovarFun>(covariance_function),
      std::forward<CovarArgs>(covar_args), ops, msgs);
}

}  // namespace math
}  // namespace stan

#endif
