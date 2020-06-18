#ifndef STAN_MATH_PRIM_PROB_POISSON_BINOMIAL_LPMF_HPP
#define STAN_MATH_PRIM_PROB_POISSON_BINOMIAL_LPMF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/poisson_binomial_log_probs.hpp>


namespace stan {
namespace math {

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
template<bool propto, typename T_y, typename T_theta>
return_type_t<T_theta> poisson_binomial_lpmf(
  const T_y &y, const T_theta &theta) {
  static const char *function = "poisson_binomial_lpmf";

  int size_theta = size_mvt(theta);
  if (size_theta > 1) {
    check_consistent_sizes(function, "Successes variables", y,
                           "Probability parameters", theta);
  }

  scalar_seq_view<T_y> y_vec(y);
  vector_seq_view < T_theta > theta_vec(theta);

  int max_sz = max_size(y, size_theta);
  std::cout << max_sz << "\n";
  for (int i = 0; i < max_sz; ++i) {
    check_bounded(function, "Successes variable", y_vec[i], 0, theta_vec[i].size());
    check_finite(function, "Probability parameters", theta_vec[i]);
    check_bounded(function, "Probability parameters", theta_vec[i], 0.0, 1.0);
    std::cout << "i\n" << y_vec[i] << "\n" << theta_vec[i] << "\n\n";
  }

  return_type_t<T_theta> log_prob = 0.0;
  for (int i = 0; i < max_sz; ++i) {
    auto x = poisson_binomial_log_probs(y_vec[i], theta_vec[i]);
    log_prob += x(y_vec[i]);
  }
  return log_prob;
}

template<typename T_y, typename T_theta>
return_type_t<T_theta> poisson_binomial_lpmf(
  const T_y &y, const T_theta &theta) {
  return poisson_binomial_lpmf<false>(y, theta);
}

}  // namespace math
}  // namespace stan
#endif
