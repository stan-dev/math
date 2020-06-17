#ifndef STAN_MATH_PRIM_PROB_POISSON_BINOMIAL_LCDF_HPP
#define STAN_MATH_PRIM_PROB_POISSON_BINOMIAL_LCDF_HPP

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
#include <stan/math/prim/fun/poisson_binomial_log_probs.hpp>


namespace stan {
namespace math {

/** \ingroup prob_dists
 * Returns the log CDF for the Poisson-binomial distribution evaluated at the
 * specified number of successes and probabilities of successes.
 *
 * @tparam T_theta type of chance of success parameters
 * @param y array of numbers of successes
 * @param theta array of chances of success parameters
 * @return sum of log probabilities
 * @throw std::domain_error if y is out of bounds
 * @throw std::domain_error if theta is not a valid vector of probabilities
 * @throw std::invalid_argument If y and theta are different lengths
 */
template <bool propto, typename T_theta>
return_type_t<T_theta> poisson_binomial_lcdf(
    const std::vector<int>& y, const std::vector<Eigen::Matrix<T_theta, Eigen::Dynamic, 1> >& theta) {
  static const char* function = "poisson_binomial_lcdf";

  check_consistent_sizes(function, "Successes variables", y,
                         "Probability parameters", theta);

  for (int i = 0; i < theta.size(); ++i) {
    check_bounded(function, "Successes variable", y[i], 0, theta[i].size());
    check_finite(function, "Probability parameters", theta[i]);
    check_bounded(function, "Probability parameters", theta[i], 0.0, 1.0);
  }

  T_theta P = 0.0;
  for (int i = 0; i < theta.size(); ++i) {
    auto alpha = poisson_binomial_log_probs(y[i], theta[i]);
    std::vector<T_theta> Pi(y[i] + 1);
    for (int j = 0; j <= y[i]; ++j) {
      Pi[j] = alpha(j);
    }
    P += log_sum_exp(Pi);
  }

  return P;
}

template <typename T_theta>
return_type_t<T_theta> poisson_binomial_lcdf(
    const std::vector<int>& y, const std::vector<Eigen::Matrix<T_theta, Eigen::Dynamic, 1> >& theta) {
  return poisson_binomial_lcdf<false>(y, theta);
}

/** \ingroup prob_dists
 * Returns the log CDF for the Poisson-binomial distribution evaluated at the
 * specified number of successes and probabilities of successes.
 *
 * @tparam T_theta type of chance of success parameters
 * @param y array of numbers of successes
 * @param theta chance of success parameters
 * @return sum of log probabilities
 * @throw std::domain_error if y is out of bounds
 * @throw std::domain_error if theta is not a valid vector of probabilities
 */
template <bool propto, typename T_theta>
return_type_t<T_theta> poisson_binomial_lcdf(const std::vector<int>& y,
                                             const Eigen::Matrix<T_theta, Eigen::Dynamic, 1>& theta) {
  static const char* function = "poisson_binomial_lcdf";

  check_bounded(function, "Successes variable", y, 0, theta.size());
  check_finite(function, "Probability parameters", theta);
  check_bounded(function, "Probability parameters", theta, 0.0, 1.0);

  // find largest y and build log-prob matrix only for largest y
  int max_y = *std::max_element(y.begin(), y.end());
  auto alpha = poisson_binomial_log_probs(max_y, theta);

  T_theta P = 0.0;
  for (int i = 0; i < y.size(); ++i) {
    std::vector<T_theta> Pi(y[i] + 1);
    for (int j = 0; j <= y[i]; ++j) {
      Pi[j] = alpha(j);
    }
    P += log_sum_exp(Pi);
  }

  return P;
}

template <typename T_theta>
inline return_type_t<T_theta> poisson_binomial_lcdf(const std::vector<int>& y,
                                                    const Eigen::Matrix<T_theta, Eigen::Dynamic, 1>& theta) {
  return poisson_binomial_lcdf<false>(y, theta);
}

/** \ingroup prob_dists
 * Returns the log CDF for the Poisson-binomial distribution evaluated at the
 * specified number of successes and probabilities of successes.
 *
 * @tparam T_theta type of chance of success parameters
 * @param y number of successes
 * @param theta chance of success parameters
 * @return log probability
 * @throw std::domain_error if y is out of bounds
 * @throw std::domain_error if theta is not a valid vector of probabilities
 */
template <bool propto, typename T_theta>
return_type_t<T_theta> poisson_binomial_lcdf(int y, const Eigen::Matrix<T_theta, Eigen::Dynamic, 1>& theta) {
  const std::vector<int> ys{y};
  return poisson_binomial_lcdf<propto>(ys, theta);
}

template <typename T_theta>
inline return_type_t<T_theta> poisson_binomial_lcdf(int y,
                                                    const Eigen::Matrix<T_theta, Eigen::Dynamic, 1>& theta) {
  return poisson_binomial_lcdf<false>(y, theta);
}

}  // namespace math
}  // namespace stan
#endif  // STAN_MATH_PRIM_PROB_POISSON_BINOMIAL_LCDF_HPP
