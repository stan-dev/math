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
 * @tparam ThetaVec A type inheriting from `Eigen::EigenBase` with dynamic
 * sized rows and 1 column.
 * @tparam CovarFun A functor with an
 *  `operator()(CovarArgsElements...)` method. The `operator()` method should accept as
 * arguments the inner elements of `CovarArgs`. The return type of the
 * `operator()` method should be a type inheriting from `Eigen::EigenBase` with
 * dynamic sized rows and columns.
 * @tparam CovarArgs A tuple of types to passed as the first arguments of
 * `CovarFun::operator()`
 * @param[in] y observed counts
 * @param[in] y_index group to which each observation belongs
 * @param[in] theta_0 the initial guess for the Laplace approximation.
 * @param covariance_function Function that returns covariance function.
 * @param covar_args arguments for the covariance function.
 * @param tolerance Tolerated gradient norm for Newton solver.
 * @param max_num_steps maximum number of steps before the Newton solver
 *                      breaks and returns an error.
 * @param hessian_block_size the size of the block for a block-diagonal
 *              Hessian of the log likelihood, i.e second derivative of
 *              log p(y|theta,phi) wrt theta.
 * @param solver Type of Newton solver. Each corresponds to a distinct choice
 *               of B matrix (i.e. application SWM formula):
 *               1. computes square-root of negative Hessian.
 *               2. computes square-root of covariance matrix.
 *               3. computes no square-root and uses LU decomposition.
 * @param max_steps_line_search Number of steps after which the algorithm
 *                        gives up on doing a linesearch. If 0, no linesearch.
 * @param[in] max_steps_line_search
 * @param msgs message stream for the covariance and likelihood function.
 */
template <bool propto = false, typename ThetaVec, typename CovarFun, typename CovarArgs,
          require_all_eigen_vector_t<ThetaVec>* = nullptr>
inline auto laplace_marginal_tol_poisson_log_lpmf(
    const std::vector<int>& y, const std::vector<int>& y_index,
    const ThetaVec& theta_0, CovarFun&& covariance_function,
    CovarArgs&& covar_args, double tolerance, int64_t max_num_steps,
    const int hessian_block_size, const int solver,
    const int max_steps_line_search, std::ostream* msgs) {
  laplace_options ops{hessian_block_size, solver, max_steps_line_search,
                      tolerance, max_num_steps};
  return laplace_marginal_density(
      poisson_log_likelihood{}, std::forward_as_tuple(y, y_index), theta_0,
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
 * @tparam ThetaVec A type inheriting from `Eigen::EigenBase` with dynamic
 * sized rows and 1 column.
 * @tparam CovarFun A functor with an
 *  `operator()(CovarArgsElements...)` method. The `operator()` method should accept as
 * arguments the inner elements of `CovarArgs`. The return type of the
 * `operator()` method should be a type inheriting from `Eigen::EigenBase` with
 * dynamic sized rows and columns.
 * @tparam CovarArgs A tuple of types to passed as the first arguments of
 * `CovarFun::operator()`
 * @param[in] y observed counts
 * @param[in] y_index group to which each observation belongs
 * @param[in] theta_0 the initial guess for the Laplace approximation.
 * @param covariance_function Function that returns covariance function.
 * @param covar_args arguments for the covariance function.
 * @param msgs message stream for the covariance and likelihood function.
 */
template <bool propto = false, typename ThetaVec, typename CovarFun, typename CovarArgs,
          require_eigen_vector_t<ThetaVec>* = nullptr>
inline auto laplace_marginal_poisson_log_lpmf(const std::vector<int>& y,
                                              const std::vector<int>& y_index,
                                              const ThetaVec& theta_0,
                                              CovarFun&& covariance_function,
                                              CovarArgs&& covar_args,
                                              std::ostream* msgs) {
  constexpr laplace_options ops{1, 1, 0, 1e-6, 100};
  return laplace_marginal_density(
      poisson_log_likelihood{}, std::forward_as_tuple(y, y_index), theta_0,
      covariance_function, std::forward<CovarArgs>(covar_args), ops, msgs);
}

}  // namespace math
}  // namespace stan

#endif
