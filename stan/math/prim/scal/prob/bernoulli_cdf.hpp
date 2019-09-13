#ifndef STAN_MATH_PRIM_SCAL_PROB_BERNOULLI_CDF_HPP
#define STAN_MATH_PRIM_SCAL_PROB_BERNOULLI_CDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_bounded.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/fun/size_zero.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Returns the CDF of the Bernoulli distribution. If containers are
 * supplied, returns the product of the probabilities.
 *
 * @tparam T_n type of integer parameter
 * @tparam T_prob type of chance of success parameter
 * @param n integer parameter
 * @param theta chance of success parameter
 * @return probability or product of probabilities
 * @throw std::domain_error if theta is not a valid probability
 * @throw std::invalid_argument if container sizes mismatch.
 */
template <typename T_n, typename T_prob>
inline auto bernoulli_cdf(T_n&& n, T_prob&& theta) {
  using T_partials = partials_return_t<T_n, T_prob>;
  T_partials P(1.0);

  static const char* function = "bernoulli_cdf";
  check_finite(function, "Probability parameter", theta);
  check_bounded(function, "Probability parameter", theta, 0.0, 1.0);
  check_consistent_sizes(function, "Random variable", n,
                         "Probability parameter", theta);

  const scalar_seq_view<T_n> n_vec(n);
  const scalar_seq_view<T_prob> theta_vec(theta);
  const size_t size = max_size(n, theta);
  operands_and_partials<T_prob> ops_partials(theta);
  if (size_zero(n, theta)) {
    return ops_partials.build(P);
  }

  for (size_t i = 0; i < size; i++) {
    const auto n_val = value_of(n_vec[i]);
    // Explicit return for extreme values
    // The gradients are technically ill-defined, but treated as zero
    if (n_val < 0) {
      return ops_partials.build(T_partials(1.0));
    }
    // Explicit results for extreme values
    // The gradients are technically ill-defined, but treated as zero
    if (n_val >= 1) {
      continue;
    }
    const T_partials Pi = 1 - value_of(theta_vec[i]);
    P *= Pi;
    if (!is_constant_all<T_prob>::value) {
      ops_partials.edge1_.partials_[i] += -1 / Pi;
    }
  }
  if (!is_constant_all<T_prob>::value) {
    for (size_t i = 0; i < stan::length(theta); ++i) {
      ops_partials.edge1_.partials_[i] *= P;
    }
  }
  return ops_partials.build(P);
}

}  // namespace math
}  // namespace stan
#endif
