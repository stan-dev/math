#ifndef STAN_MATH_LAPLACE_LAPLACE_APPROX_BERNOULLI_RNG_HPP
#define STAN_MATH_LAPLACE_LAPLACE_APPROX_BERNOULLI_RNG_HPP

#include <stan/math/mix/laplace/prob/laplace_base_rng.hpp>

namespace stan {
namespace math {

/**
 * In a latent gaussian model,
 *
 *   theta ~ Normal(theta | 0, Sigma(phi))
 *   y ~ pi(y | theta)
 *
 * return a multivariate normal random variate sampled
 * from the gaussian approximation of p(theta | y, phi),
 * where the likelihood is a Bernoulli with logit link.
 */
template <typename CovarFun, typename ThetaMatrix, class RNG,
          typename TrainTuple, typename PredTuple, typename... Args,
          require_eigen_t<ThetaMatrix>* = nullptr>
inline Eigen::VectorXd  // CHECK -- right return type
laplace_marginal_tol_bernoulli_logit_rng(
    const std::vector<int>& y, const std::vector<int>& n_samples,
    const double tolerance, const long int max_num_steps,
    const int hessian_block_size, const int solver,
    const int max_steps_line_search, const ThetaMatrix& theta_0,
    CovarFun&& covariance_function, RNG& rng, std::ostream* msgs,
    TrainTuple&& train_tuple, PredTuple&& pred_tuple, Args&&... args) {
  Eigen::VectorXd eta_dummy(0);
  return laplace_base_rng(
      diff_bernoulli_logit(to_vector(n_samples), to_vector(y)),
      covariance_function, eta_dummy, theta_0, rng, msgs, tolerance,
      max_num_steps, hessian_block_size, solver, max_steps_line_search,
      std::forward<TrainTuple>(train_tuple),
      std::forward<PredTuple>(pred_tuple), std::forward<Args>(args)...);
}

template <typename CovarFun, typename ThetaMatrix, class RNG,
          typename TrainTuple, typename PredTuple, typename... Args,
          require_eigen_t<ThetaMatrix>* = nullptr>
inline Eigen::VectorXd  // CHECK -- right return type
laplace_marginal_bernoulli_logit_rng(const std::vector<int>& y,
                                     const std::vector<int>& n_samples,
                                     const ThetaMatrix& theta_0,
                                     CovarFun&& covariance_function, RNG& rng,
                                     std::ostream* msgs,
                                     TrainTuple&& train_tuple,
                                     PredTuple&& pred_tuple, Args&&... args) {
  constexpr double tolerance = 1e-6;
  constexpr long int max_num_steps = 100;
  constexpr int hessian_block_size = 1;
  constexpr int solver = 1;
  constexpr int max_steps_line_search = 0;
  Eigen::VectorXd eta_dummy(0);
  return laplace_base_rng(
      diff_bernoulli_logit(to_vector(n_samples), to_vector(y)),
      covariance_function, eta_dummy, theta_0, rng, msgs, tolerance,
      max_num_steps, hessian_block_size, solver, max_steps_line_search,
      std::forward<TrainTuple>(train_tuple),
      std::forward<PredTuple>(pred_tuple), std::forward<Args>(args)...);
}

}  // namespace math
}  // namespace stan

#endif
