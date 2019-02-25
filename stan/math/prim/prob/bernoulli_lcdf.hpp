#ifndef STAN_MATH_PRIM_PROB_BERNOULLI_LCDF_HPP
#define STAN_MATH_PRIM_PROB_BERNOULLI_LCDF_HPP

#include <stanh/prim/meta/is_constant_struct.hpp>
#include <stanh/prim/meta/partials_return_type.hpp>
#include <stanh/prim/meta/operands_and_partials.hpp>
#include <stanh/prim/err/check_consistent_sizes.hpp>
#include <stanh/prim/err/check_bounded.hpp>
#include <stanh/prim/err/check_finite.hpp>
#include <stanh/prim/err/check_not_nan.hpp>
#include <stanh/prim/fun/size_zero.hpp>
#include <stanh/prim/fun/constants.hpp>
#include <stanh/prim/fun/inv_logit.hpp>
#include <stanh/prim/fun/log1m.hpp>
#include <stanh/prim/fun/value_of.hpp>
#include <stanh/prim/meta/include_summand.hpp>
#include <stanh/prim/metaar_seq_view.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Returns the log CDF of the Bernoulli distribution. If containers are
 * supplied, returns the log sum of the probabilities.
 *
 * @tparam T_n type of integer parameter
 * @tparam T_prob type of chance of success parameter
 * @param n integer parameter
 * @param theta chance of success parameter
 * @return log probability or log sum of probabilities
 * @throw std::domain_error if theta is not a valid probability
 * @throw std::invalid_argument if container sizes mismatch.
 */
template <typename T_n, typename T_prob>
typename return_type<T_prob>::type bernoulli_lcdf(const T_n& n,
                                                  const T_prob& theta) {
  static const char* function = "bernoulli_lcdf";
  typedef
      typename stan::partials_return_type<T_n, T_prob>::type T_partials_return;

  if (size_zero(n, theta))
    return 0.0;

  T_partials_return P(0.0);

  check_finite(function, "Probability parameter", theta);
  check_bounded(function, "Probability parameter", theta, 0.0, 1.0);
  check_consistent_sizes(function, "Random variable", n,
                         "Probability parameter", theta);

  scalar_seq_view<T_n> n_vec(n);
  scalar_seq_view<T_prob> theta_vec(theta);
  size_t size = max_size(n, theta);

  using std::log;
  operands_and_partials<T_prob> ops_partials(theta);

  // Explicit return for extreme values
  // The gradients are technically ill-defined, but treated as zero
  for (size_t i = 0; i < stan::length(n); i++) {
    if (value_of(n_vec[i]) < 0)
      return ops_partials.build(negative_infinity());
  }

  for (size_t i = 0; i < size; i++) {
    // Explicit results for extreme values
    // The gradients are technically ill-defined, but treated as zero
    if (value_of(n_vec[i]) >= 1)
      continue;

    const T_partials_return Pi = 1 - value_of(theta_vec[i]);

    P += log(Pi);

    if (!is_constant_struct<T_prob>::value)
      ops_partials.edge1_.partials_[i] -= 1 / Pi;
  }

  return ops_partials.build(P);
}

}  // namespace math
}  // namespace stan
#endif
