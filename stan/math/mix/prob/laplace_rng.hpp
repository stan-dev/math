#ifndef STAN_MATH_MIX_PROB_LAPLACE_RNG_HPP
#define STAN_MATH_MIX_PROB_LAPLACE_RNG_HPP

#include <stan/math/mix/functor/laplace_base_rng.hpp>
#include <stan/math/mix/functor/laplace_likelihood.hpp>

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
 * where the log likelihood is given by L_f.
 */
template <typename LFun, typename EtaVec, typename DeltaVec, typename CovarFun,
          typename ThetaVec, typename RNG, typename TrainTuple,
          typename PredTuple, typename... Args>
inline Eigen::VectorXd laplace_marginal_tol_rng(
    LFun&& L_f, const EtaVec& eta, const DeltaVec& delta_L,
    const std::vector<int>& delta_int_L, const double tolerance,
    const long int max_num_steps, const int hessian_block_size,
    const int solver, const int max_steps_line_search, const ThetaVec& theta_0,
    CovarFun&& K_f, RNG& rng, std::ostream* msgs, TrainTuple&& train_tuple,
    PredTuple&& pred_tuple, Args&&... args) {
  laplace_options ops{hessian_block_size, solver, max_steps_line_search,
                      tolerance, max_num_steps};
  return laplace_base_rng(std::forward<LFun>(L_f),
      std::forward_as_tuple(delta_L, delta_int_L),
      K_f, eta, theta_0, rng, msgs, ops, std::forward<TrainTuple>(train_tuple),
      std::forward<PredTuple>(pred_tuple), std::forward<Args>(args)...);
}

template <typename LFun, typename DeltaVec, typename CovarFun,
          typename ThetaVec, typename RNG, typename TrainTuple,
          typename PredTuple, typename... Args>
inline Eigen::VectorXd laplace_marginal_tol_rng(
    LFun&& L_f, const DeltaVec& delta_L, const std::vector<int>& delta_int_L,
    const double tolerance, const long int max_num_steps,
    const int hessian_block_size, const int solver,
    const int max_steps_line_search, const ThetaVec& theta_0, CovarFun&& K_f,
    RNG& rng, std::ostream* msgs, TrainTuple&& train_tuple,
    PredTuple&& pred_tuple, Args&&... args) {
  laplace_options ops{hessian_block_size, solver, max_steps_line_search,
                      tolerance, max_num_steps};
  Eigen::Matrix<double, 0, 0> eta;
  return laplace_base_rng(std::forward<LFun>(L_f),
      std::forward_as_tuple(delta_L, delta_int_L),
      K_f, eta, theta_0, rng, msgs, ops, std::forward<TrainTuple>(train_tuple),
      std::forward<PredTuple>(pred_tuple), std::forward<Args>(args)...);
}

template <typename LFun, typename EtaVec, typename DeltaVec, typename CovarFun,
          typename ThetaVec, typename RNG, typename TrainTuple,
          typename PredTuple, typename... Args>
inline Eigen::VectorXd laplace_marginal_rng(
    LFun&& L_f, const EtaVec& eta, const DeltaVec& delta_L,
    const std::vector<int>& delta_int_L, const ThetaVec& theta_0,
    CovarFun&& K_f, RNG& rng, std::ostream* msgs, TrainTuple&& train_tuple,
    PredTuple&& pred_tuple, Args&&... args) {
  constexpr laplace_options ops{1, 1, 0, 1e-6, 100};
  return laplace_base_rng(std::forward<LFun>(L_f),
      std::forward_as_tuple(delta_L, delta_int_L),
      K_f, eta, theta_0, rng, msgs, ops, std::forward<TrainTuple>(train_tuple),
      std::forward<PredTuple>(pred_tuple), std::forward<Args>(args)...);
}

template <typename LFun, typename DeltaVec, typename CovarFun,
          typename ThetaVec, typename RNG, typename TrainTuple,
          typename PredTuple, typename... Args>
inline Eigen::VectorXd laplace_marginal_rng(
    LFun&& L_f, const DeltaVec& delta_L, const std::vector<int>& delta_int_L,
    const ThetaVec& theta_0, CovarFun&& K_f, RNG& rng, std::ostream* msgs,
    TrainTuple&& train_tuple, PredTuple&& pred_tuple, Args&&... args) {
  constexpr laplace_options ops{1, 1, 0, 1e-6, 100};
  Eigen::Matrix<double, 0, 0> eta;
  return laplace_base_rng(std::forward<LFun>(L_f),
      std::forward_as_tuple(delta_L, delta_int_L),
      K_f, eta, theta_0, rng, msgs, ops, std::forward<TrainTuple>(train_tuple),
      std::forward<PredTuple>(pred_tuple), std::forward<Args>(args)...);
}

}  // namespace math
}  // namespace stan

#endif
