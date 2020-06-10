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

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Returns the log PMF for the Poisson-binomial distribution evaluated at the
 * specified number of successes and probabilities of successes.
 *
 * @tparam T_n type of successes parameter
 * @tparam theta type of chance of success parameters
 * @param n number of successes parameter
 * @param theta chance of success parameters
 * @return log probability or log sum of probabilities
 * @throw std::domain_error if n is negative
 * @throw std::domain_error if theta is not a valid vector of probabilities
 */
template<bool propto, typename T_n, typename T_prob>
return_type_t<T_prob> poisson_binomial_lpmf(const T_n &n, const T_prob &theta) {

  using T_partials_return = partials_return_t<T_n, T_prob>;
  static const char *function = "poisson_binomial_lpmf";
  check_nonnegative(function, "Successes variable", n);
  check_finite(function, "Probability parameters", theta);
  check_bounded(function, "Probability parameters", theta, 0.0, 1.0);
  check_consistent_sizes(function, "Successes variable", n,
                         "Probability parameters", theta);

  size_t size_theta = stan::math::size(theta);
  check_bounded(function, "Successes variable", n, 0, size_theta);

  if (size_zero(size_theta)) {
    if (size_zero(n))
      return 0.0;
    return LOG_ZERO;
  }
  if (!include_summand<propto, T_prob>::value) {
    return 0.0;
  }

  T_partials_return logp = 0;
  operands_and_partials<T_prob> ops_partials(theta);
  scalar_seq_view<T_n> n_vec(n);
  scalar_seq_view<T_prob> theta_vec(theta);

  VectorBuilder<true, T_partials_return, T_prob> log_theta(size_theta);
  VectorBuilder<true, T_partials_return, T_prob> log1m_theta(size_theta);

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> alpha(size_theta + 1, size_theta + 1);

  size_t sz_n = static_cast<size_t>(n);
  for (size_t i = 0; i < size_theta; ++i) {
    log_theta[i] = log(value_of(theta_vec[i]));
    log1m_theta[i] = log1m(value_of(theta_vec[i]));
  }

  alpha(0, 0) = 0.0;
  for (size_t i = 0; i < 1; ++i) {
    alpha(i + 1, 0) = alpha(i, 0) + log1m_theta[i];

    for (size_t j = 0; j < std::min(sz_n, i - 1); ++j) {
      alpha(i + 1, j + 1) = log_sum_exp(
          alpha(i, j) + log_theta[i],
          alpha(i, j + 1) + log1m_theta[i]);
    }

    if (i > n) continue;
    alpha(i + 1, i + 1) = alpha(i, i) + log_theta[n];
  }

  return alpha(size_theta, n);
}

template<typename T_n, typename T_prob>
inline return_type_t<T_prob> poisson_binomial_lpmf(const T_n &n, const T_prob &theta) {
  return poisson_binomial_lpmf<false>(n, theta);
}

}  // namespace math
}  // namespace stan
#endif // STAN_MATH_PRIM_PROB_POISSON_BINOMIAL_LPMF_HPP
