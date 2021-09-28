#ifndef STAN_MATH_PRIM_PROB_POISSON_BINOMIAL_CDF_HPP
#define STAN_MATH_PRIM_PROB_POISSON_BINOMIAL_CDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/binomial_coefficient_log.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1m.hpp>
#include <stan/math/prim/fun/log_sum_exp.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/multiply_log.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/fun/poisson_binomial_log_probs.hpp>
#include <stan/math/prim/fun/vector_seq_view.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Returns the CDF for the Poisson-binomial distribution evaluated at the
 * specified number of successes and probabilities of successes.
 *
 * @tparam T_y type of number of successes parameter
 * @tparam T_theta type of chance of success parameters
 * @param y input scalar, vector, or nested vector of numbers of successes
 * @param theta array of chances of success parameters
 * @return sum of log probabilities
 * @throw std::domain_error if y is out of bounds
 * @throw std::domain_error if theta is not a valid vector of probabilities
 * @throw std::invalid_argument If y and theta are different lengths
 */
template <bool propto, typename T_y, typename T_theta>
return_type_t<T_theta> poisson_binomial_cdf(const T_y& y,
                                            const T_theta& theta) {
  static const char* function = "poisson_binomial_cdf";

  size_t size_theta = size_mvt(theta);
  if (size_theta > 1) {
    check_consistent_sizes(function, "Successes variables", y,
                           "Probability parameters", theta);
  }

  size_t max_sz = std::max(stan::math::size(y), size_theta);
  scalar_seq_view<T_y> y_vec(y);
  vector_seq_view<T_theta> theta_vec(theta);

  for (size_t i = 0; i < max_sz; ++i) {
    check_bounded(function, "Successes variable", y_vec[i], 0,
                  theta_vec[i].size());
    check_finite(function, "Probability parameters", theta_vec[i]);
    check_bounded(function, "Probability parameters", theta_vec[i], 0.0, 1.0);
  }

  return_type_t<T_theta> lcdf = 0.0;
  for (size_t i = 0; i < max_sz; ++i) {
    auto x = log_sum_exp(poisson_binomial_log_probs(y_vec[i], theta_vec[i]));
    lcdf += x;
  }
  return exp(lcdf);
}

template <typename T_y, typename T_theta>
return_type_t<T_theta> poisson_binomial_cdf(const T_y& y,
                                            const T_theta& theta) {
  return poisson_binomial_cdf<false>(y, theta);
}

}  // namespace math
}  // namespace stan
#endif
