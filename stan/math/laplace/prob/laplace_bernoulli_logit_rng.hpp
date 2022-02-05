#ifndef STAN_MATH_LAPLACE_LAPLACE_APPROX_BERNOULLI_RNG_HPP
#define STAN_MATH_LAPLACE_LAPLACE_APPROX_BERNOULLI_RNG_HPP

#include <stan/math/laplace/prob/laplace_base_rng.hpp>

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
template <typename CovarFun, typename T_theta, typename T_x, class RNG,
          typename TupleData, typename... Args>
inline Eigen::VectorXd  // CHECK -- right return type
laplace_bernoulli_logit_rng(
    const std::vector<int>& y, const std::vector<int>& n_samples,
    CovarFun&& covariance_function,
    const Eigen::Matrix<T_theta, Eigen::Dynamic, 1>& theta_0,
    const TupleData& data_tuple, RNG& rng, std::ostream* msgs = nullptr,
    const double tolerance = 1e-6, const long int max_num_steps = 100,
    const int hessian_block_size = 0, const int solver = 1,
    const int do_line_search = 0, const int max_steps_line_search = 10,
    Args&&... args) {
  Eigen::VectorXd eta_dummy(0);
  return laplace_base_rng(
      diff_bernoulli_logit(to_vector(n_samples), to_vector(y)),
      covariance_function, eta_dummy, theta_0, data_tuple, rng, msgs, tolerance,
      max_num_steps, hessian_block_size, solver, do_line_search,
      max_steps_line_search, std::forward<Args>(args)...);
}

}  // namespace math
}  // namespace stan

#endif
