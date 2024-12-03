#ifndef STAN_MATH_MIX_FUNCTOR_LAPLACE_BASE_RNG_HPP
#define STAN_MATH_MIX_FUNCTOR_LAPLACE_BASE_RNG_HPP

#include <stan/math/mix/functor/laplace_marginal_density.hpp>
#include <stan/math/prim/fun.hpp>
#include <stan/math/prim/prob/multi_normal_cholesky_rng.hpp>

#include <Eigen/Sparse>

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
 * @tparam D Type of likelihood function.
 * @tparam LLArgs Type of arguments of likelihood function.
 * @tparam ThetaMatrix Type of latent Gaussian variables.
 * @tparam EtaMatrix Type of additional arguments for likelihood function.
 * @tparam CovarFun Type of covariance function.
 * @tparam RNG Type of RNG number.
 * @tparam TrainTuple Type of observed/training covariate.
 * @tparam PredTuple Type of predictive covariate.
 * @tparam CovarArgs Tuple of additional arguments for the covariance function
 * @param ll_fun Likelihood function.
 * @param ll_args Arguments for likelihood function.
 * @param covariance_function Covariance function.
 * @param eta Additional arguments for likelihood function.
 * @param theta_0 Initial guess for finding the mode of the conditional
                  pi(theta_pred | y, phi, x_pred).
 * @param options Control parameter for optimizer underlying Laplace approx.
 * @param train_tuple Observed/training covariates for covariance function.
 * @param pred_tuple Predictive covariates for covariance function.
 * @param rng Rng number.
 * @param msgs Stream for function prints.
 * @param args Variadic arguments for likelihood function.
 */
template <typename D, typename LLArgs, typename ThetaMatrix,
          typename CovarFun, class RNG, typename TrainTuple, typename PredTuple,
          typename CovarArgs,
          require_all_eigen_t<ThetaMatrix>* = nullptr>
inline Eigen::VectorXd laplace_base_rng(
    D&& ll_fun, LLArgs&& ll_args, CovarFun&& covariance_function,
    const ThetaMatrix& theta_0,
    const laplace_options& options, TrainTuple&& train_tuple,
    PredTuple&& pred_tuple, RNG& rng, std::ostream* msgs,
    CovarArgs&& covar_args) {
  using Eigen::MatrixXd;
  using Eigen::VectorXd;
  using EtaMatrix = Eigen::Matrix<double, 0, 0>;
  EtaMatrix eta{};
  auto covar_args_val = stan::math::to_ref(value_of(covar_args));
  auto eta_dbl = value_of(eta);
  auto md_est = laplace_marginal_density_est(
      ll_fun, ll_args, eta_dbl, value_of(theta_0), covariance_function,
      std::tuple_cat(std::forward<TrainTuple>(train_tuple), covar_args_val),
      options, msgs);
  // Modified R&W method
  MatrixXd covariance_pred = apply(
      [&covariance_function, &msgs](auto&&... args_val) {
        return covariance_function(args_val..., msgs);
      },
      std::tuple_cat(std::forward<PredTuple>(pred_tuple), covar_args_val));
  VectorXd pred_mean = covariance_pred * md_est.theta_grad;
  if (options.solver == 1 || options.solver == 2) {
    Eigen::MatrixXd V_dec = mdivide_left_tri<Eigen::Lower>(
        md_est.L, md_est.W_r * covariance_pred);
    Eigen::MatrixXd Sigma = covariance_pred - V_dec.transpose() * V_dec;
    return multi_normal_rng(pred_mean, Sigma, rng);
  } else {
    Eigen::MatrixXd Sigma
        = covariance_pred
          - covariance_pred
                * (md_est.W_r
                   - md_est.W_r
                         * md_est.LU.solve(md_est.covariance * md_est.W_r))
                * covariance_pred;
    return multi_normal_rng(pred_mean, Sigma, rng);
  }
}

}  // namespace math
}  // namespace stan

#endif
