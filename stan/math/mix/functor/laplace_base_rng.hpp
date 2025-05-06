#ifndef STAN_MATH_MIX_FUNCTOR_LAPLACE_BASE_RNG_HPP
#define STAN_MATH_MIX_FUNCTOR_LAPLACE_BASE_RNG_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/mix/functor/laplace_marginal_density.hpp>
#include <stan/math/prim/prob/multi_normal_cholesky_rng.hpp>



namespace stan {
namespace math {

/**
 * In a latent gaussian model,
 *
 *   theta ~ Normal(theta | 0, Sigma(phi, x))
 *   y ~ pi(y | theta, eta)
 *
 * returns a multivariate normal random variate sampled
 * from the Laplace approximation of p(theta_pred | y, phi, x_pred).
 * Note that while the data is observed at x (train_tuple), the new samples
 * are drawn for covariates x_pred (pred_tuple).
 * To sample the "original" theta's, set pred_tuple = train_tuple.
 * @tparam LLFunc Type of likelihood function.
 * @tparam LLArgs Tuple of arguments types of likelihood function.
 * @tparam ThetaVec A type inheriting from `Eigen::EigenBase` with dynamic
 * sized rows and 1 column.
 * @tparam CovarFun A functor with an
 *  `operator()(CovarArgsElements..., {TrainTupleElements...|
 PredTupleElements...})`
 *  method. The `operator()` method should accept as arguments the
 *  inner elements of `CovarArgs`. The return type of the `operator()` method
 *  should be a type inheriting from `Eigen::EigenBase` with dynamic sized
 *  rows and columns.
 * @tparam RNG A valid boost rng type
 * @tparam CovarArgs A tuple of types to passed as the first arguments of
 `CovarFun::operator()`
 * @param ll_fun Likelihood function.
 * @param ll_args Arguments for likelihood function.
 * @param theta_0 Initial guess for finding the mode of the conditional
                  pi(theta_pred | y, phi, x_pred).
 * @param covariance_function Covariance function.
 * @param covar_args Observed/training covariates for covariance function.
 * @param options Control parameter for optimizer underlying Laplace approx.
 * @param rng Rng number.
 * @param msgs Stream for function prints.
 */
template <
    typename LLFunc, typename LLArgs, typename ThetaVec, typename CovarFun,
    typename CovarArgs, typename RNG, require_all_eigen_t<ThetaVec>* = nullptr,
    require_t<is_all_arithmetic_scalar<CovarArgs, LLArgs, ThetaVec>>* = nullptr>
inline Eigen::VectorXd laplace_base_rng(LLFunc&& ll_fun, LLArgs&& ll_args,
                                        ThetaVec&& theta_0,
                                        CovarFun&& covariance_function,
                                        CovarArgs&& covar_args,
                                        const laplace_options& options,
                                        RNG& rng, std::ostream* msgs) {
  auto md_est = laplace_marginal_density_est(
      ll_fun, std::forward<LLArgs>(ll_args), std::forward<ThetaVec>(theta_0),
      std::forward<CovarFun>(covariance_function),
      to_ref(std::forward<CovarArgs>(covar_args)), options, msgs);
  // Modified R&W method
  auto&& covariance_train = md_est.covariance;
  Eigen::VectorXd mean_train = covariance_train * md_est.theta_grad;
  if (options.solver == 1 || options.solver == 2) {
    Eigen::MatrixXd V_dec
        = md_est.L.template triangularView<Eigen::Lower>().solve(
            md_est.W_r * covariance_train);
    Eigen::MatrixXd Sigma = covariance_train - V_dec.transpose() * V_dec;
    return multi_normal_rng(std::move(mean_train), std::move(Sigma), rng);
  } else {
    Eigen::MatrixXd Sigma
        = covariance_train
          - covariance_train
                * (md_est.W_r
                   - md_est.W_r
                         * md_est.LU.solve(covariance_train * md_est.W_r))
                * covariance_train;
    return multi_normal_rng(std::move(mean_train), std::move(Sigma), rng);
  }
}

}  // namespace math
}  // namespace stan

#endif
