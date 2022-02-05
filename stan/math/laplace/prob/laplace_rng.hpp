#ifndef STAN_MATH_LAPLACE_LAPLACE_RNG_HPP
#define STAN_MATH_LAPLACE_LAPLACE_RNG_HPP

#include <stan/math/laplace/prob/laplace_rng.hpp>

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
template <typename LFun, typename T_eta, typename CovarFun, typename T_theta,
          typename RNG, typename TupleData, typename... Args>
inline Eigen::VectorXd laplace_rng(
    LFun&& L_f, const Eigen::Matrix<T_eta, Eigen::Dynamic, 1>& eta,
    const Eigen::VectorXd& delta_L, const std::vector<int>& delta_int_L,
    CovarFun&& K_f, const Eigen::Matrix<T_theta, Eigen::Dynamic, 1>& theta_0,
    const TupleData& data_tuple, RNG& rng, std::ostream* msgs = nullptr,
    const double tolerance = 1e-6, const long int max_num_steps = 100,
    const int hessian_block_size = 0, const int solver = 1,
    const int do_line_search = 0, const int max_steps_line_search = 10,
    Args&&... args) {
  return laplace_base_rng(diff_likelihood<LFun>(std::forward<LFun>(L_f),
                                                delta_L, delta_int_L, msgs),
                          K_f, eta, theta_0, data_tuple, rng, msgs, tolerance,
                          max_num_steps, hessian_block_size, solver,
                          do_line_search, max_steps_line_search,
                          std::forward<Args>(args)...);
}

}  // namespace math
}  // namespace stan

#endif
