#ifndef STAN_MATH_PRIM_PROB_POISSON_BINOMIAL_LPMF_HPP
#define STAN_MATH_PRIM_PROB_POISSON_BINOMIAL_LPMF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/binomial_coefficient_log.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1m.hpp>
#include <stan/math/prim/fun/log_sum_exp.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/multiply_log.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>

#include <iostream>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Returns the log PMF for the Poisson-binomial distribution evaluated at the
 * specified number of successes and probabilities of successes.
 *
 * @tparam T_theta type of chance of success parameters
 * @param y number of successes parameter
 * @param theta chance of success parameters
 * @return log probability or log sum of probabilities
 * @throw std::domain_error if y is negative
 * @throw std::domain_error if theta is not a valid vector of probabilities
 */
template<bool propto, typename T_theta>
return_type_t<T_theta> poisson_binomial_lpmf(
  int y,
  const Eigen::Matrix<T_theta, Eigen::Dynamic, 1>& theta) {

  static const char *function = "poisson_binomial_lpmf";
  check_nonnegative(function, "Successes variable", y);
  check_bounded(function, "Successes variable", y, 0, theta.size());
  check_finite(function, "Probability parameters", theta);
  check_bounded(function, "Probability parameters", theta, 0.0, 1.0);

  size_t size_theta = theta.size();
  Eigen::Matrix<T_theta, Eigen::Dynamic, 1> log_theta(size_theta);
  Eigen::Matrix<T_theta, Eigen::Dynamic, 1> log1m_theta(size_theta);
  Eigen::Matrix<T_theta, Eigen::Dynamic, Eigen::Dynamic> alpha(
    size_theta + 1, size_theta + 1);

  for (size_t i = 0; i < size_theta; ++i) {
    log_theta[i] = log(theta(i));
    log1m_theta[i] = log1m(theta(i));
  }

  size_t sz_y = static_cast<size_t>(y);
  // alpha[i, j] = log prob of j successes in first i trials
  alpha(0, 0) = 0.0;
  for (size_t i = 0; i < size_theta; ++i) {
    // no success in i trials
    alpha(i + 1, 0) = alpha(i, 0) + log1m_theta[i];

    // 0 < j < i successes in i trials
    for (size_t j = 0; j < std::min(sz_y, i); ++j) {
      alpha(i + 1, j + 1) = log_sum_exp(
          alpha(i, j) + log_theta[i],
          alpha(i, j + 1) + log1m_theta[i]);
    }

    // i successes in i trials
    if (i < sz_y) {
      alpha(i + 1, i + 1) = alpha(i, i) + log_theta(i);
    }
  }

  std::cout << alpha << std::endl;
  return alpha(size_theta, y);
}

template<typename T_theta>
inline return_type_t<T_theta> poisson_binomial_lpmf(
  int y,
  const Eigen::Matrix<T_theta, Eigen::Dynamic, 1>& theta) {
  return poisson_binomial_lpmf<false>(y, theta);
}

}  // namespace math
}  // namespace stan
#endif // STAN_MATH_PRIM_PROB_POISSON_BINOMIAL_LPMF_HPP
