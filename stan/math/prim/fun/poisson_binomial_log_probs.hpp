#ifndef STAN_MATH_PRIM_FUN_POISSON_BINOMIAL_LOG_PROBS_HPP
#define STAN_MATH_PRIM_FUN_POISSON_BINOMIAL_LOG_PROBS_HPP

#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1m.hpp>
#include <stan/math/prim/fun/log_sum_exp.hpp>


namespace stan {
namespace math {

/**
 * Return a row vector of ones
 *
 * @param K size of the row vector
 * @return A row vector of size K with all elements initialised to 1.
 * @throw std::domain_error if K is negative.
 */
template <typename T_theta>
Eigen::Matrix<T_theta, Eigen::Dynamic, Eigen::Dynamic> poisson_binomial_log_probs(
  int y, const Eigen::Matrix<T_theta, Eigen::Dynamic, 1>& theta) {

  int size_theta = theta.size();
  using vec = Eigen::Matrix<T_theta, Eigen::Dynamic, 1>;

  vec log_theta(size_theta);
  vec log1m_theta(size_theta);

  for (int i = 0; i < size_theta; ++i) {
    log_theta[i] = log(theta(i));
    log1m_theta[i] = log1m(theta(i));
  }

  Eigen::Matrix<T_theta, Eigen::Dynamic, Eigen::Dynamic> alpha(size_theta + 1, y + 1);

  // alpha[i, j] = log prob of j successes in first i trials
  alpha(0, 0) = 0.0;
  for (int i = 0; i < size_theta; ++i) {
    // no success in i trials
    alpha(i + 1, 0) = alpha(i, 0) + log1m_theta[i];

    // 0 < j < i successes in i trials
    for (int j = 0; j < std::min(y, i); ++j) {
      alpha(i + 1, j + 1) = log_sum_exp(alpha(i, j) + log_theta[i],
                                        alpha(i, j + 1) + log1m_theta[i]);
    }

    // i successes in i trials
    if (i < y) {
      alpha(i + 1, i + 1) = alpha(i, i) + log_theta(i);
    }
  }

  return alpha.row(size_theta);
}

}  // namespace math
}  // namespace stan

#endif // STAN_MATH_PRIM_FUN_POISSON_BINOMIAL_LOG_PROBS_HPP
