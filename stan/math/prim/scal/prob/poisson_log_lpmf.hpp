#ifndef STAN_MATH_PRIM_SCAL_PROB_POISSON_LOG_LPMF_HPP
#define STAN_MATH_PRIM_SCAL_PROB_POISSON_LOG_LPMF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_nonnegative.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/fun/size_zero.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/fun/lgamma.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <cmath>
#include <limits>

namespace stan {
namespace math {

// PoissonLog(n|alpha)  [n >= 0]   = Poisson(n|exp(alpha))
template <bool propto, typename T_n, typename T_log_rate>
inline auto poisson_log_lpmf(T_n&& n, T_log_rate&& alpha) {
  using T_partials = partials_return_t<T_n, T_log_rate>;
  T_partials logp(0.0);
  using std::exp;
  static const char* function = "poisson_log_lpmf";
  check_nonnegative(function, "Random variable", n);
  check_not_nan(function, "Log rate parameter", alpha);
  check_consistent_sizes(function, "Random variable", n, "Log rate parameter",
                         alpha);

  const scalar_seq_view<T_n> n_vec(n);
  const scalar_seq_view<T_log_rate> alpha_vec(alpha);
  const size_t size = max_size(n, alpha);
  operands_and_partials<T_log_rate> ops_partials(alpha);
  if (!include_summand<propto, T_log_rate>::value) {
    return ops_partials.build(logp);
  } else if (size_zero(n, alpha)) {
    return ops_partials.build(logp);
  }

  for (size_t i = 0; i < size; i++) {
    if (std::numeric_limits<double>::infinity() == alpha_vec[i]) {
      return ops_partials.build(T_partials(LOG_ZERO));
    }
    if (-std::numeric_limits<double>::infinity() == alpha_vec[i]
        && n_vec[i] != 0) {
      return ops_partials.build(T_partials(LOG_ZERO));
    }
    const auto alpha_val = value_of(alpha_vec[i]);
    const auto n_val = n_vec[i];
    T_partials exp_alpha = 0;
    if (include_summand<propto>::value) {
      logp -= lgamma(n_val + 1.0);
    }
    if (include_summand<propto, T_log_rate>::value) {
      exp_alpha = exp(alpha_val);
      logp += n_val * alpha_val - exp_alpha;
    }
    if (!is_constant_all<T_log_rate>::value) {
      ops_partials.edge1_.partials_[i] += n_val - exp_alpha;
    }
  }
  return ops_partials.build(logp);
}

template <typename T_n, typename T_log_rate>
inline auto poisson_log_lpmf(T_n&& n, T_log_rate&& alpha) {
  return poisson_log_lpmf<false>(n, alpha);
}

}  // namespace math
}  // namespace stan
#endif
