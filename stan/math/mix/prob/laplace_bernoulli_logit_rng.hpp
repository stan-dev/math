#ifndef STAN_MATH_MIX_PROB_LAPLACE_BERNOULLI_LOGIT_RNG_HPP
#define STAN_MATH_MIX_PROB_LAPLACE_BERNOULLI_LOGIT_RNG_HPP

#include <stan/math/mix/functor/laplace_base_rng.hpp>
#include <stan/math/mix/functor/laplace_likelihood.hpp>
#include <stan/math/mix/prob/laplace_marginal_bernoulli_logit_lpmf.hpp>

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
    const ThetaMatrix& theta_0, CovarFun&& covariance_function,
    TrainTuple&& train_tuple, PredTuple&& pred_tuple, const double tolerance,
    const int64_t max_num_steps, const int hessian_block_size,
    const int solver, const int max_steps_line_search, RNG& rng,
    std::ostream* msgs, Args&&... args) {
  laplace_options ops{hessian_block_size, solver, max_steps_line_search,
                      tolerance, max_num_steps};
  Eigen::Matrix<double, 0, 0> eta_dummy;
  return laplace_base_rng(
      bernoulli_logit_likelihood{},
      std::forward_as_tuple(to_vector(y), n_samples), covariance_function,
      eta_dummy, theta_0, rng, msgs, ops, std::forward<TrainTuple>(train_tuple),
      std::forward<PredTuple>(pred_tuple), std::forward<Args>(args)...);
}

template <typename CovarFun, typename ThetaMatrix, class RNG,
          typename TrainTuple, typename PredTuple, typename... Args,
          require_eigen_t<ThetaMatrix>* = nullptr>
inline Eigen::VectorXd  // CHECK -- right return type
laplace_marginal_bernoulli_logit_rng(const std::vector<int>& y,
                                     const std::vector<int>& n_samples,
                                     const ThetaMatrix& theta_0,
                                     CovarFun&& covariance_function,
                                     TrainTuple&& train_tuple,
                                     PredTuple&& pred_tuple, RNG& rng,
                                     std::ostream* msgs, Args&&... args) {
  constexpr laplace_options ops{1, 1, 0, 1e-6, 100};
  Eigen::Matrix<double, 0, 0> eta_dummy;
  return laplace_base_rng(bernoulli_logit_likelihood{},
                          std::forward_as_tuple(to_vector(y), n_samples),
                          covariance_function, eta_dummy, theta_0, ops,
                          std::make_tuple(), std::make_tuple(), rng, msgs,
                          std::forward<Args>(args)...);
}

}  // namespace math
}  // namespace stan

#endif
