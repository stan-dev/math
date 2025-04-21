#ifndef STAN_MATH_PRIM_FUN_POISSON_BINOMIAL_LOG_PROBS_HPP
#define STAN_MATH_PRIM_FUN_POISSON_BINOMIAL_LOG_PROBS_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1m.hpp>
#include <stan/math/prim/fun/log_sum_exp.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/vector_seq_view.hpp>

namespace stan {
namespace math {

/**
 * Returns the last row of the log probability matrix of the Poisson-Binomial
 * distribution given the number of successes and a vector of success
 * probabilities.
 *
 * @tparam T_theta template expression
 * @param y numbers of successes
 * @param theta N-dimensional vector of success probabilities for each trial
 * @return the last row of the computed log probability matrix
 */
template <typename T_theta, typename T_scalar = scalar_type_t<T_theta>,
          require_vector_t<T_theta>* = nullptr>
Eigen::Matrix<T_scalar, 1, -1>
  poisson_binomial_log_probs(int y, const T_theta& theta) {
  int size_theta = theta.size();
  plain_type_t<T_theta> log_theta = log(theta);
  plain_type_t<T_theta> log1m_theta = log1m(theta);

  Eigen::Matrix<T_scalar, Eigen::Dynamic, Eigen::Dynamic> alpha(size_theta + 1,
                                                                y + 1);

  // alpha[i, j] = log prob of j successes in first i trials
  alpha(0, 0) = 0.0;
  for (int i = 0; i < size_theta; ++i) {
    // no success in i trials
    alpha(i + 1, 0) = alpha(i, 0) + log1m_theta[i];

    // 0 < j < i successes in i trials
    for (int j = 0; j < std::min(y, i); ++j) {
      alpha(i + 1, j + 1) = log_mix(theta[i], alpha(i, j), alpha(i, j + 1));
    }

    // i successes in i trials
    if (i < y) {
      alpha(i + 1, i + 1) = alpha(i, i) + log_theta[i];
    }
  }

  return alpha.row(size_theta);
}

}  // namespace math
}  // namespace stan
#endif
