#ifndef STAN_MATH_PRIM_PROB_POISSON_LCCDF_HPP
#define STAN_MATH_PRIM_PROB_POISSON_LCCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/log_gamma_q_dgamma.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>
#include <limits>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * The log complementary cumulative distribution function of the Poisson
 * distribution, log P(n + 1, lambda) with the regularized lower incomplete
 * gamma function P. It is computed in log space, so it stays finite where
 * the probability underflows.
 *
 * @tparam T_n type of n
 * @tparam T_rate type of lambda
 * @param n random variable
 * @param lambda rate parameter
 * @return log probability, or sum of log probabilities for containers
 * @throw std::domain_error if lambda is negative
 * @throw std::invalid_argument if container sizes mismatch
 */
template <typename T_n, typename T_rate>
inline return_type_t<T_rate> poisson_lccdf(const T_n& n, const T_rate& lambda) {
  using T_partials_return = partials_return_t<T_n, T_rate>;
  using T_n_ref = ref_type_if_not_constant_t<T_n>;
  using T_lambda_ref = ref_type_if_not_constant_t<T_rate>;
  static constexpr const char* function = "poisson_lccdf";
  check_consistent_sizes(function, "Random variable", n, "Rate parameter",
                         lambda);

  T_n_ref n_ref = n;
  T_lambda_ref lambda_ref = lambda;

  decltype(auto) lambda_val
      = to_ref(as_value_column_array_or_scalar(lambda_ref));

  check_nonnegative(function, "Rate parameter", lambda_val);

  if (size_zero(n, lambda)) {
    return 0;
  }

  auto ops_partials = make_partials_propagator(lambda_ref);

  scalar_seq_view<T_n_ref> n_vec(n_ref);
  scalar_seq_view<T_lambda_ref> lambda_vec(lambda_ref);
  const size_t N = max_size(n, lambda);

  // An element with probability 0 (n >= 0 with lambda = 0, or n = INT_MAX,
  // which stands for an infinite count) makes the result -inf; the
  // gradients are treated as zero
  for (size_t i = 0; i < N; i++) {
    if (n_vec.val(i) >= 0
        && (value_of_rec(lambda_vec.val(i)) == 0
            || n_vec.val(i) == std::numeric_limits<int>::max())) {
      return ops_partials.build(negative_infinity());
    }
  }

  T_partials_return log_P(0.0);
  for (size_t i = 0; i < N; i++) {
    // An element with probability 1 (n < 0) contributes nothing
    if (n_vec.val(i) < 0) {
      continue;
    }
    const T_partials_return lambda_dbl = lambda_vec.val(i);
    const auto tail
        = internal::log_gamma_tail<false, false, is_autodiff_v<T_rate>>(
            n_vec.val(i) + 1.0, lambda_dbl);
    log_P += tail.log_prob;
    if constexpr (is_autodiff_v<T_rate>) {
      partials<0>(ops_partials)[i] += tail.d_z;
    }
  }
  return ops_partials.build(log_P);
}

}  // namespace math
}  // namespace stan
#endif
