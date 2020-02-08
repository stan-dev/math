#ifndef STAN_MATH_PRIM_PROB_POISSON_LOG_LPMF_HPP
#define STAN_MATH_PRIM_PROB_POISSON_LOG_LPMF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <cmath>

namespace stan {
namespace math {

// PoissonLog(n|alpha)  [n >= 0]   = Poisson(n|exp(alpha))
template <bool propto, typename T_n, typename T_log_rate>
return_type_t<T_log_rate> poisson_log_lpmf(const T_n& n,
                                           const T_log_rate& alpha) {
  using T_partials_return = partials_return_t<T_n, T_log_rate>;

  static const char* function = "poisson_log_lpmf";

  using std::exp;

  if (size_zero(n, alpha)) {
    return 0.0;
  }

  T_partials_return logp(0.0);

  check_nonnegative(function, "Random variable", n);
  check_not_nan(function, "Log rate parameter", alpha);
  check_consistent_sizes(function, "Random variable", n, "Log rate parameter",
                         alpha);

  if (!include_summand<propto, T_log_rate>::value) {
    return 0.0;
  }

  scalar_seq_view<T_n> n_vec(n);
  scalar_seq_view<T_log_rate> alpha_vec(alpha);
  size_t max_size_seq_view = max_size(n, alpha);

  for (size_t i = 0, size_alpha = size(alpha); i < size_alpha; i++) {
    if (INFTY == alpha_vec[i]) {
      return LOG_ZERO;
    }
  }
  for (size_t i = 0; i < max_size_seq_view; i++) {
    if (NEGATIVE_INFTY == alpha_vec[i] && n_vec[i] != 0) {
      return LOG_ZERO;
    }
  }

  operands_and_partials<T_log_rate> ops_partials(alpha);

  VectorBuilder<include_summand<propto>::value, T_partials_return, T_n>
      lgamma_n_plus_one(size(n));
  if (include_summand<propto>::value) {
    for (size_t i = 0, size_n = size(n); i < size_n; i++) {
      lgamma_n_plus_one[i] = lgamma(n_vec[i] + 1.0);
    }
  }

  VectorBuilder<include_summand<propto, T_log_rate>::value, T_partials_return,
                T_log_rate>
      exp_alpha(size(alpha));
  for (size_t i = 0, size_alpha = size(alpha); i < size_alpha; i++) {
    exp_alpha[i] = exp(value_of(alpha_vec[i]));
  }

  for (size_t i = 0; i < max_size_seq_view; i++) {
    const auto& alpha_val = value_of(alpha_vec[i]);
    if (!(alpha_val == NEGATIVE_INFTY && n_vec[i] == 0)) {
      if (include_summand<propto>::value) {
        logp -= lgamma_n_plus_one[i];
      }
      logp += n_vec[i] * alpha_val - exp_alpha[i];
    }

    if (!is_constant_all<T_log_rate>::value) {
      ops_partials.edge1_.partials_[i] += n_vec[i] - exp_alpha[i];
    }
  }
  return ops_partials.build(logp);
}

template <typename T_n, typename T_log_rate>
inline return_type_t<T_log_rate> poisson_log_lpmf(const T_n& n,
                                                  const T_log_rate& alpha) {
  return poisson_log_lpmf<false>(n, alpha);
}

}  // namespace math
}  // namespace stan
#endif
