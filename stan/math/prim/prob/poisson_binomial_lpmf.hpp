#ifndef STAN_MATH_PRIM_PROB_POISSON_BINOMIAL_LPMF_HPP
#define STAN_MATH_PRIM_PROB_POISSON_BINOMIAL_LPMF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1m.hpp>
#include <stan/math/prim/fun/log_sum_exp.hpp>
#include <stan/math/prim/fun/size.hpp>

namespace stan {
namespace math {

template<typename T>
using vec = Eigen::Matrix<T, Eigen::Dynamic, 1>;
template<typename T>
using mat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

template<typename T_theta>
mat<T_theta> compute_alpha_(int y, const vec<T_theta>& theta) {

  size_t sz_y = static_cast<size_t>(y);
  size_t size_theta = theta.size();
  vec<T_theta> log_theta(size_theta);
  vec<T_theta> log1m_theta(size_theta);

  for (size_t i = 0; i < size_theta; ++i) {
    log_theta[i] = log(theta(i));
    log1m_theta[i] = log1m(theta(i));
  }

  mat<T_theta> alpha(size_theta + 1, sz_y + 1);

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

  return alpha;
}

/** \ingroup prob_dists
 * Returns the log PMF for the Poisson-binomial distribution evaluated at an
 * specified array of numbers of successes and probabilities of successes.
 *
 * @tparam T_theta type of chance of success parameters
 * @param y array of numbers of successes
 * @param theta array of chances of success parameters
 * @return sum of log probabilities
 * @throw std::domain_error if y is out of bounds
 * @throw std::domain_error if theta is not a valid vector of probabilities
 * @throw std::invalid_argument If y and theta are different lengths
 */
template<bool propto, typename T_theta>
return_type_t<T_theta> poisson_binomial_lpmf(
  const std::vector<int>& y,
  const std::vector< vec<T_theta> >& theta) {
  static const char *function = "poisson_binomial_lpmf";

  int sz_theta = static_cast<int>(theta.size());
  check_consistent_sizes(function,
    "Successes variables", y,
    "Probability parameters", theta);

  for (int i = 0; i < sz_theta; ++i) {
    check_bounded(function, "Successes variable", y[i], 0, theta[i].size());
    check_finite(function, "Probability parameters", theta[i]);
    check_bounded(function, "Probability parameters", theta[i], 0.0, 1.0);
  }

  T_theta log_prob = 0.0;
  for (int i = 0; i < sz_theta; ++i) {
    mat<T_theta> alpha = compute_alpha_(y[i], theta[i]);
    log_prob += alpha(theta[i].size(), y[i]);
  }

  return log_prob;
}

template<typename T_theta>
return_type_t<T_theta> poisson_binomial_lpmf(
  const std::vector<int>& y,
  const std::vector< vec<T_theta> >& theta) {
  return poisson_binomial_lpmf<false>(y, theta);
}

/** \ingroup prob_dists
 * Returns the log PMF for the Poisson-binomial distribution evaluated at an
 * specified array of numbers of successes and probabilities of successes.
 *
 * @tparam T_theta type of chance of success parameters
 * @param y array of numbers of successes
 * @param theta chance of success parameters
 * @return sum of log probabilities
 * @throw std::domain_error if y is out of bounds
 * @throw std::domain_error if theta is not a valid vector of probabilities
 */
template<bool propto, typename T_theta>
return_type_t<T_theta> poisson_binomial_lpmf(
  const std::vector<int>& y,
  const vec<T_theta>& theta) {
  static const char *function = "poisson_binomial_lpmf";

  check_bounded(function, "Successes variable", y, 0, theta.size());
  check_finite(function, "Probability parameters", theta);
  check_bounded(function, "Probability parameters", theta, 0.0, 1.0);

  // find largest y and build log-prob matrix only for largest y
  int max_y = *std::max_element(y.begin(), y.end());
  mat<T_theta> alpha = compute_alpha_(max_y, theta);

  T_theta log_prob = 0;
  for (std::vector<int>::size_type i = 0; i < y.size(); ++i) {
    log_prob += alpha(theta.size(), y[i]);
  }

  return log_prob;
}

template<typename T_theta>
inline return_type_t<T_theta> poisson_binomial_lpmf(
  const std::vector<int>& y,
  const vec<T_theta>& theta) {
  return poisson_binomial_lpmf<false>(y, theta);
}

/** \ingroup prob_dists
 * Returns the log PMF for the Poisson-binomial distribution evaluated at the
 * specified number of successes and probabilities of successes.
 *
 * @tparam T_theta type of chance of success parameters
 * @param y number of successes
 * @param theta chance of success parameters
 * @return log probability
 * @throw std::domain_error if y is out of bounds
 * @throw std::domain_error if theta is not a valid vector of probabilities
 */
template<bool propto, typename T_theta>
return_type_t<T_theta> poisson_binomial_lpmf(
  int y,
  const vec<T_theta>& theta) {
  const std::vector<int> ys{y};
  return poisson_binomial_lpmf<propto>(ys, theta);
}

template<typename T_theta>
inline return_type_t<T_theta> poisson_binomial_lpmf(
  int y,
  const vec<T_theta>& theta) {
  return poisson_binomial_lpmf<false>(y, theta);
}

}  // namespace math
}  // namespace stan
#endif // STAN_MATH_PRIM_PROB_POISSON_BINOMIAL_LPMF_HPP
