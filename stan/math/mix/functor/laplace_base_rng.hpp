#ifndef STAN_MATH_MIX_FUNCTOR_LAPLACE_BASE_RNG_HPP
#define STAN_MATH_MIX_FUNCTOR_LAPLACE_BASE_RNG_HPP

#include <stan/math/mix/functor/laplace_marginal_density.hpp>
#include <stan/math/prim/prob/multi_normal_cholesky_rng.hpp>
#include <stan/math/prim/fun.hpp>

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
 * from the gaussian approximation of p(theta_pred | y, phi, x_pred).
 * Note that while the data is observed at x, the new samples
 * are drawn for covariates x_pred.
 * To sample the "original" theta's, set x_pred = x.
 */
template <typename D, typename LLArgs, typename ThetaMatrix, typename EtaMatrix,
          typename CovarFun, class RNG, typename TrainTuple, typename PredTuple,
          typename... Args,
          require_all_eigen_t<ThetaMatrix, EtaMatrix>* = nullptr>
inline Eigen::VectorXd laplace_base_rng(
    D&& ll_fun, LLArgs&& ll_args, CovarFun&& covariance_function,
    const ThetaMatrix& eta, const EtaMatrix& theta_0,
    const laplace_options& options, TrainTuple&& train_tuple,
    PredTuple&& pred_tuple, RNG& rng, std::ostream* msgs, Args&&... args) {
  using Eigen::MatrixXd;
  using Eigen::VectorXd;
  auto args_dbl = std::make_tuple(to_ref(value_of(args))...);
  auto eta_dbl = value_of(eta);
  auto md_est = apply(
      [&](auto&&... args_val) {
        return laplace_marginal_density_est(
            ll_fun, ll_args, covariance_function, eta_dbl, value_of(theta_0),
            msgs, options, args_val...);
      },
      std::tuple_cat(std::forward<TrainTuple>(train_tuple), args_dbl));
  // Modified R&W method
  MatrixXd covariance_pred = apply(
      [&covariance_function, &msgs](auto&&... args_val) {
        return covariance_function(args_val..., msgs);
      },
      std::tuple_cat(std::forward<PredTuple>(pred_tuple), args_dbl));
  VectorXd pred_mean = covariance_pred * md_est.l_grad.head(theta_0.rows());
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
