#ifndef STAN_MATH_PRIM_PROB_POISSON_LPMF_HPP
#define STAN_MATH_PRIM_PROB_POISSON_LPMF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/fun/is_inf.hpp>
#include <stan/math/prim/scal/fun/lgamma.hpp>
#include <stan/math/prim/scal/fun/multiply_log.hpp>
#include <stan/math/prim/scal/fun/size_zero.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>

namespace stan {
namespace math {

// Poisson(n|lambda)  [lambda > 0;  n >= 0]
template <bool propto, typename T_n, typename T_rate>
return_type_t<T_rate> poisson_lpmf(const T_n& n, const T_rate& lambda) {
  using T_partials_return = partials_return_t<T_n, T_rate>;

  static const char* function = "poisson_lpmf";

  if (size_zero(n, lambda)) {
    return 0.0;
  }

  T_partials_return logp(0.0);

  check_nonnegative(function, "Random variable", n);
  check_not_nan(function, "Rate parameter", lambda);
  check_nonnegative(function, "Rate parameter", lambda);
  check_consistent_sizes(function, "Random variable", n, "Rate parameter",
                         lambda);

  if (!include_summand<propto, T_rate>::value) {
    return 0.0;
  }

  scalar_seq_view<T_n> n_vec(n);
  scalar_seq_view<T_rate> lambda_vec(lambda);
  size_t max_size_seq_view = max_size(n, lambda);

  for (size_t i = 0, size_lambda = size(lambda); i < size_lambda; i++) {
    if (is_inf(lambda_vec[i])) {
      return LOG_ZERO;
    }
  }
  for (size_t i = 0; i < max_size_seq_view; i++) {
    if (lambda_vec[i] == 0 && n_vec[i] != 0) {
      return LOG_ZERO;
    }
  }

  operands_and_partials<T_rate> ops_partials(lambda);

  VectorBuilder<include_summand<propto>::value, T_partials_return, T_n>
      lgamma_n_plus_one(size(n));
  if (include_summand<propto>::value) {
    for (size_t i = 0, size_n = size(n); i < size_n; i++) {
      lgamma_n_plus_one[i] = lgamma(n_vec[i] + 1.0);
    }
  }

  for (size_t i = 0; i < max_size_seq_view; i++) {
    const auto& lambda_val = value_of(lambda_vec[i]);
    if (!(lambda_val == 0 && n_vec[i] == 0)) {
      if (include_summand<propto>::value) {
        logp -= lgamma_n_plus_one[i];
      }
      logp += multiply_log(n_vec[i], lambda_val) - lambda_val;
    }

    if (!is_constant_all<T_rate>::value) {
      ops_partials.edge1_.partials_[i] += n_vec[i] / lambda_val - 1.0;
    }
  }

  return ops_partials.build(logp);
}

template <typename T_n, typename T_rate>
inline return_type_t<T_rate> poisson_lpmf(const T_n& n, const T_rate& lambda) {
  return poisson_lpmf<false>(n, lambda);
}

}  // namespace math
}  // namespace stan
#endif
