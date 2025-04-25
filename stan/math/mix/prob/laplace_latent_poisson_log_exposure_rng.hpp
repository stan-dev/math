#ifndef STAN_MATH_MIX_PROB_LAPLACE_LATENT_POISSON_LOG_EXPOSURE_RNG_HPP
#define STAN_MATH_MIX_PROB_LAPLACE_LATENT_POISSON_LOG_EXPOSURE_RNG_HPP

#include <stan/math/mix/functor/laplace_base_rng.hpp>
#include <stan/math/mix/functor/laplace_likelihood.hpp>
#include <stan/math/mix/prob/laplace_marginal_poisson_log_exposure_lpmf.hpp>

namespace stan {
namespace math {

/**
 * Wrapper function around the laplace_marginal function for
 * a log poisson likelihood with exposure. Returns the marginal density
 * p(y | phi) by marginalizing out the latent gaussian variable,
 * with a Laplace approximation. See the laplace_marginal function
 * for more details.
 *
 * @tparam ThetaVec A type inheriting from `Eigen::EigenBase` with dynamic
 * sized rows and 1 column.
 * @tparam CovarFun A functor with an
 *  `operator()(CovarArgsElements...)` method. The `operator()` method should
 * accept as arguments the inner elements of `CovarArgs`. The return type of the
 * `operator()` method should be a type inheriting from `Eigen::EigenBase` with
 * dynamic sized rows and columns.
 * @tparam CovarArgs A tuple of types to passed as the first arguments of
 * `CovarFun::operator()`
 * @tparam RNG A valid boost rng type
 * @param y
 * @param y_index
 * @param ye
 * @param theta_0
 * @param covariance_function
 * @param covar_args
 * @param tolerance
 * @param max_num_steps
 * @param hessian_block_size
 * @param solver
 * @param max_steps_line_search
 * @param rng
 * @param msgs
 *
 */
template <typename YeVec, typename ThetaVec, typename CovarFun,
          typename CovarArgs, typename RNG,
          require_eigen_t<ThetaVec>* = nullptr>
inline auto  // CHECK -- right return type
laplace_latent_tol_poisson_2_log_rng(
    const std::vector<int>& y, const std::vector<int>& y_index, const YeVec& ye,
    ThetaVec&& theta_0, CovarFun&& covariance_function, CovarArgs&& covar_args,
    const double tolerance, const int64_t max_num_steps,
    const int hessian_block_size, const int solver,
    const int max_steps_line_search, RNG& rng, std::ostream* msgs) {
  laplace_options ops{hessian_block_size, solver, max_steps_line_search,
                      tolerance, max_num_steps};
  return laplace_base_rng(poisson_log_exposure_likelihood{},
                          std::forward_as_tuple(y, y_index, ye),
                          std::forward<ThetaVec>(theta_0),
                          std::forward<CovarFun>(covariance_function),
                          std::forward<CovarArgs>(covar_args), ops, rng, msgs);
}

/**
 * Wrapper function around the laplace_marginal function for
 * a log poisson likelihood with exposure. Returns the marginal density
 * p(y | phi) by marginalizing out the latent gaussian variable,
 * with a Laplace approximation. See the laplace_marginal function
 * for more details.
 *
 * @tparam ThetaVec A type inheriting from `Eigen::EigenBase` with dynamic
 * sized rows and 1 column.
 * @tparam CovarFun A functor with an
 *  `operator()(CovarArgsElements...)` method. The `operator()` method should
 * accept as arguments the inner elements of `CovarArgs`. The return type of the
 * `operator()` method should be a type inheriting from `Eigen::EigenBase` with
 * dynamic sized rows and columns.
 * @tparam CovarArgs A tuple of types to passed as the first arguments of
 * `CovarFun::operator()`
 * @tparam RNG A valid boost rng type
 * @param y
 * @param y_index
 * @param ye
 * @param theta_0
 * @param covariance_function
 * @param covar_args
 * @param rng
 * @param msgs
 *
 */
template <typename YeVec, typename ThetaVec, typename CovarFun,
          typename CovarArgs, typename RNG,
          require_eigen_t<ThetaVec>* = nullptr>
inline auto  // TODO(Steve): Allow scalar or std vector return
laplace_latent_poisson_2_log_rng(const std::vector<int>& y,
                                 const std::vector<int>& y_index,
                                 const YeVec& ye, ThetaVec&& theta_0,
                                 CovarFun&& covariance_function,
                                 CovarArgs&& covar_args, RNG& rng,
                                 std::ostream* msgs) {
  constexpr laplace_options ops{1, 1, 0, 1e-6, 100};
  return laplace_base_rng(poisson_log_exposure_likelihood{},
                          std::forward_as_tuple(y, y_index, ye),
                          std::forward<ThetaVec>(theta_0),
                          std::forward<CovarFun>(covariance_function),
                          std::forward<CovarArgs>(covar_args), ops, rng, msgs);
}

}  // namespace math
}  // namespace stan

#endif
