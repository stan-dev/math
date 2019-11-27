#ifndef STAN_MATH_PRIM_SCAL_PROB_PARETO_TYPE_2_LPDF_HPP
#define STAN_MATH_PRIM_SCAL_PROB_PARETO_TYPE_2_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_greater_or_equal.hpp>
#include <stan/math/prim/scal/err/check_nonnegative.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <stan/math/prim/scal/fun/size_zero.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <stan/math/prim/scal/fun/log1p.hpp>
#include <cmath>

namespace stan {
namespace math {

// pareto_type_2(y|lambda, alpha)  [y >= 0;  lambda > 0;  alpha > 0]
template <bool propto, typename T_y, typename T_loc, typename T_scale,
          typename T_shape>
return_type_t<T_y, T_loc, T_scale, T_shape> pareto_type_2_lpdf(
    const T_y& y, const T_loc& mu, const T_scale& lambda,
    const T_shape& alpha) {
  static const char* function = "pareto_type_2_lpdf";
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale, T_shape>;

  using std::log;

  if (size_zero(y, mu, lambda, alpha)) {
    return 0.0;
  }

  T_partials_return logp(0.0);

  check_greater_or_equal(function, "Random variable", y, mu);
  check_not_nan(function, "Random variable", y);
  check_positive_finite(function, "Scale parameter", lambda);
  check_positive_finite(function, "Shape parameter", alpha);
  check_consistent_sizes(function, "Random variable", y, "Scale parameter",
                         lambda, "Shape parameter", alpha);

  if (!include_summand<propto, T_y, T_loc, T_scale, T_shape>::value) {
    return 0.0;
  }

  scalar_seq_view<T_y> y_vec(y);
  scalar_seq_view<T_loc> mu_vec(mu);
  scalar_seq_view<T_scale> lambda_vec(lambda);
  scalar_seq_view<T_shape> alpha_vec(alpha);
  size_t N = max_size(y, mu, lambda, alpha);

  operands_and_partials<T_y, T_loc, T_scale, T_shape> ops_partials(
      y, mu, lambda, alpha);

  VectorBuilder<include_summand<propto, T_y, T_loc, T_scale, T_shape>::value,
                T_partials_return, T_y, T_loc, T_scale>
      log1p_scaled_diff(N);
  if (include_summand<propto, T_y, T_loc, T_scale, T_shape>::value) {
    for (size_t n = 0; n < N; n++) {
      log1p_scaled_diff[n] = log1p((value_of(y_vec[n]) - value_of(mu_vec[n]))
                                   / value_of(lambda_vec[n]));
    }
  }

  VectorBuilder<include_summand<propto, T_scale>::value, T_partials_return,
                T_scale>
      log_lambda(length(lambda));
  if (include_summand<propto, T_scale>::value) {
    for (size_t n = 0; n < length(lambda); n++) {
      log_lambda[n] = log(value_of(lambda_vec[n]));
    }
  }

  VectorBuilder<include_summand<propto, T_shape>::value, T_partials_return,
                T_shape>
      log_alpha(length(alpha));
  if (include_summand<propto, T_shape>::value) {
    for (size_t n = 0; n < length(alpha); n++) {
      log_alpha[n] = log(value_of(alpha_vec[n]));
    }
  }

  VectorBuilder<!is_constant_all<T_shape>::value, T_partials_return, T_shape>
      inv_alpha(length(alpha));
  if (!is_constant_all<T_shape>::value) {
    for (size_t n = 0; n < length(alpha); n++) {
      inv_alpha[n] = 1 / value_of(alpha_vec[n]);
    }
  }

  for (size_t n = 0; n < N; n++) {
    const T_partials_return y_dbl = value_of(y_vec[n]);
    const T_partials_return mu_dbl = value_of(mu_vec[n]);
    const T_partials_return lambda_dbl = value_of(lambda_vec[n]);
    const T_partials_return alpha_dbl = value_of(alpha_vec[n]);
    const T_partials_return sum_dbl = lambda_dbl + y_dbl - mu_dbl;
    const T_partials_return inv_sum = 1.0 / sum_dbl;
    const T_partials_return alpha_div_sum = alpha_dbl / sum_dbl;
    const T_partials_return deriv_1_2 = inv_sum + alpha_div_sum;

    if (include_summand<propto, T_shape>::value) {
      logp += log_alpha[n];
    }
    if (include_summand<propto, T_scale>::value) {
      logp -= log_lambda[n];
    }
    if (include_summand<propto, T_y, T_scale, T_shape>::value) {
      logp -= (alpha_dbl + 1.0) * log1p_scaled_diff[n];
    }

    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_[n] -= deriv_1_2;
    }
    if (!is_constant_all<T_loc>::value) {
      ops_partials.edge2_.partials_[n] += deriv_1_2;
    }
    if (!is_constant_all<T_scale>::value) {
      ops_partials.edge3_.partials_[n]
          -= alpha_div_sum * (mu_dbl - y_dbl) / lambda_dbl + inv_sum;
    }
    if (!is_constant_all<T_shape>::value) {
      ops_partials.edge4_.partials_[n] += inv_alpha[n] - log1p_scaled_diff[n];
    }
  }
  return ops_partials.build(logp);
}

template <typename T_y, typename T_loc, typename T_scale, typename T_shape>
inline return_type_t<T_y, T_loc, T_scale, T_shape> pareto_type_2_lpdf(
    const T_y& y, const T_loc& mu, const T_scale& lambda,
    const T_shape& alpha) {
  return pareto_type_2_lpdf<false>(y, mu, lambda, alpha);
}

}  // namespace math
}  // namespace stan
#endif
