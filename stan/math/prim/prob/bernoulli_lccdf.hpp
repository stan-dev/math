#ifndef STAN_MATH_PRIM_PROB_BERNOULLI_LCCDF_HPP
#define STAN_MATH_PRIM_PROB_BERNOULLI_LCCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Returns the log CCDF of the Bernoulli distribution. If containers are
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
return_type_t<T_prob> bernoulli_lccdf(const T_n& n, const T_prob& theta) {
  static const char* function = "bernoulli_lccdf";
  using T_partials_return = partials_return_t<T_n, T_prob>;

  if (size_zero(n, theta)) {
    return 0.0;
  }

  T_partials_return P(0.0);

  check_finite(function, "Probability parameter", theta);
  check_bounded(function, "Probability parameter", theta, 0.0, 1.0);
  check_consistent_sizes(function, "Random variable", n,
                         "Probability parameter", theta);

  scalar_seq_view<T_n> n_vec(n);
  scalar_seq_view<T_prob> theta_vec(theta);
  size_t max_size_seq_view = max_size(n, theta);

  using std::log;
  operands_and_partials<T_prob> ops_partials(theta);

  // Explicit return for extreme values
  // The gradients are technically ill-defined, but treated as zero
  for (size_t i = 0; i < size(n); i++) {
    if (value_of(n_vec[i]) < 0) {
      return ops_partials.build(0.0);
    }
  }

  for (size_t i = 0; i < max_size_seq_view; i++) {
    // Explicit results for extreme values
    // The gradients are technically ill-defined, but treated as zero
    if (value_of(n_vec[i]) >= 1) {
      return ops_partials.build(negative_infinity());
    } else {
      const T_partials_return Pi = value_of(theta_vec[i]);

      P += log(Pi);

      if (!is_constant_all<T_prob>::value) {
        ops_partials.edge1_.partials_[i] += 1 / Pi;
      }
    }
  }

  return ops_partials.build(P);
}
}  // namespace math
}  // namespace stan
#endif
