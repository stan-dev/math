#ifndef STAN_MATH_LAPLACE_PROB_LAPLACE_BASE_RNG_HPP
#define STAN_MATH_LAPLACE_PROB_LAPLACE_BASE_RNG_HPP

#include <stan/math/prim/prob/multi_normal_cholesky_rng.hpp>
#include <stan/math/prim/fun/cholesky_decompose.hpp>
#include <stan/math/laplace/laplace_marginal.hpp>

#include <Eigen/Sparse>
#include <Eigen/LU>

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
template <typename ThetaMatrix, typename EtaMatrix, typename D, typename CovarFun,
          class RNG, typename TrainTuple, typename PredTuple, typename... Args,
          require_all_eigen_t<ThetaMatrix, EtaMatrix>* = nullptr>
inline Eigen::VectorXd  // CHECK -- right return type -- It's not this need to
                        // return a std::vector<> :(
laplace_base_rng(D&& diff_likelihood, CovarFun&& covariance_function,
                 const ThetaMatrix& eta,
                 const EtaMatrix& theta_0,
                 RNG& rng,
                 std::ostream* msgs = nullptr, const double tolerance = 1e-6,
                 const long int max_num_steps = 100,
                 const int hessian_block_size = 0, const int solver = 1,
                 const int max_steps_line_search = 0,
                 TrainTuple&& train_tuple = std::tuple<>(),
                PredTuple&& pred_tuple = std::tuple<>(), Args&&... args) {
  using Eigen::MatrixXd;
  using Eigen::VectorXd;

  auto args_dbl = std::make_tuple(to_ref(value_of(args))...);

  VectorXd eta_dbl = value_of(eta);
  auto marginal_density_est = apply(
      [&](auto&&... args_val) {
        return laplace_marginal_density_est(
            diff_likelihood, covariance_function, eta_dbl, value_of(theta_0),
            msgs, tolerance, max_num_steps, hessian_block_size, solver,
            max_steps_line_search, args_val...);
      },
      std::tuple_cat(std::forward<TrainTuple>(train_tuple), args_dbl));
  auto marginal_density = marginal_density_est.lmd;
  Eigen::SparseMatrix<double> W_r = std::move(marginal_density_est.W_r);
  MatrixXd L = std::move(marginal_density_est.L);
  MatrixXd K_root = std::move(marginal_density_est.K_root);
  Eigen::PartialPivLU<MatrixXd> LU = std::move(marginal_density_est.LU);
  VectorXd l_grad = std::move(marginal_density_est.l_grad);
  MatrixXd covariance = std::move(marginal_density_est.covariance);

  // Modified R&W method
  MatrixXd covariance_pred = apply(
      [&covariance_function, &msgs](auto&&... args_val) {
        return covariance_function(args_val..., msgs);
      },
      std::tuple_cat(std::forward<PredTuple>(pred_tuple), args_dbl));

  VectorXd pred_mean = covariance_pred * l_grad.head(theta_0.rows());

  Eigen::MatrixXd Sigma;
  if (solver == 1 || solver == 2) {
    Eigen::MatrixXd V_dec
        = mdivide_left_tri<Eigen::Lower>(L, W_r * covariance_pred);
    Sigma = covariance_pred - V_dec.transpose() * V_dec;
  } else {
    Sigma = covariance_pred
            - covariance_pred * (W_r - W_r * LU.solve(covariance * W_r))
                  * covariance_pred;
  }

  return multi_normal_rng(pred_mean, Sigma, rng);
}

}  // namespace math
}  // namespace stan

#endif
