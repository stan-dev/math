#ifndef STAN_MATH_MIX_PROB_LAPLACE_POISSON_LOG_RNG_HPP
#define STAN_MATH_MIX_PROB_LAPLACE_POISSON_LOG_RNG_HPP

#include <stan/math/mix/functor/laplace_base_rng.hpp>
#include <stan/math/mix/functor/laplace_likelihood.hpp>
#include <stan/math/mix/prob/laplace_marginal_poisson_log_exposure_lpmf.hpp>
#include <stan/math/mix/prob/laplace_marginal_poisson_log_lpmf.hpp>

namespace stan {
namespace math {

/**
 * In a latent gaussian model,
 *
 *   theta ~ Normal(theta | 0, Sigma(phi))
 *   y ~ pi(y | theta)
 *
 * return a multivariate normal random variate sampled
 * from the gaussian approximation of p(theta | y, phi)
 * where the likelihood is a Poisson with a log link.
 */
template <typename CovarFun, typename ThetaMatrix, class RNG,
          typename TrainTuple, typename PredTuple, typename... Args,
          require_eigen_t<ThetaMatrix>* = nullptr>
inline Eigen::VectorXd laplace_marginal_tol_poisson_log_rng(
    const std::vector<int>& y, const std::vector<int>& n_samples,
    const double tolerance, const long int max_num_steps,
    const int hessian_block_size, const int solver,
    const int max_steps_line_search, const ThetaMatrix& theta_0,
    CovarFun&& covariance_function, RNG& rng, std::ostream* msgs,
    TrainTuple&& train_tuple, PredTuple&& pred_tuple, Args&&... args) {
  Eigen::VectorXd eta_dummy;
  poisson_log_likelihood L;
  laplace_options ops{hessian_block_size, solver,
    max_steps_line_search, tolerance, max_num_steps};
  return laplace_base_rng(
      laplace_likelihood<poisson_log_likelihood>(L, to_vector(y), n_samples, msgs),
      covariance_function, eta_dummy, theta_0, rng, msgs, ops,
      std::forward<TrainTuple>(train_tuple),
      std::forward<PredTuple>(pred_tuple), std::forward<Args>(args)...);
}


template <typename CovarFun, typename ThetaMatrix, class RNG,
          typename TrainTuple, typename PredTuple, typename... Args,
          require_eigen_t<ThetaMatrix>* = nullptr>
inline Eigen::VectorXd laplace_marginal_poisson_log_rng(
    const std::vector<int>& y, const std::vector<int>& n_samples,
    const ThetaMatrix& theta_0, CovarFun&& covariance_function, RNG& rng,
    std::ostream* msgs, TrainTuple&& train_tuple, PredTuple&& pred_tuple,
    Args&&... args) {
  Eigen::VectorXd eta_dummy;
  constexpr laplace_options ops{1, 1, 0, 1e-6, 100};
  return laplace_base_rng(
      laplace_likelihood<poisson_log_likelihood>(poisson_log_likelihood{},
                                              to_vector(y), n_samples, msgs),
      covariance_function, eta_dummy, theta_0, rng, msgs, ops,
      std::forward<TrainTuple>(train_tuple),
      std::forward<PredTuple>(pred_tuple), std::forward<Args>(args)...);
}

}  // namespace math
}  // namespace stan

#endif
