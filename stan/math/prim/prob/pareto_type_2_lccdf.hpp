#ifndef STAN_MATH_PRIM_PROB_PARETO_TYPE_2_LCCDF_HPP
#define STAN_MATH_PRIM_PROB_PARETO_TYPE_2_LCCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T_y, typename T_loc, typename T_scale, typename T_shape>
return_type_t<T_y, T_loc, T_scale, T_shape> pareto_type_2_lccdf(
    const T_y& y, const T_loc& mu, const T_scale& lambda,
    const T_shape& alpha) {
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale, T_shape>;
  static const char* function = "pareto_type_2_lccdf";

  if (size_zero(y, mu, lambda, alpha)) {
    return 0.0;
  }

  using std::log;

  T_partials_return P(0.0);

  check_greater_or_equal(function, "Random variable", y, mu);
  check_not_nan(function, "Random variable", y);
  check_nonnegative(function, "Random variable", y);
  check_positive_finite(function, "Scale parameter", lambda);
  check_positive_finite(function, "Shape parameter", alpha);
  check_consistent_sizes(function, "Random variable", y, "Scale parameter",
                         lambda, "Shape parameter", alpha);

  scalar_seq_view<T_y> y_vec(y);
  scalar_seq_view<T_loc> mu_vec(mu);
  scalar_seq_view<T_scale> lambda_vec(lambda);
  scalar_seq_view<T_shape> alpha_vec(alpha);
  size_t N = max_size(y, mu, lambda, alpha);

  operands_and_partials<T_y, T_loc, T_scale, T_shape> ops_partials(
      y, mu, lambda, alpha);

  for (size_t n = 0; n < N; n++) {
    const T_partials_return y_dbl = value_of(y_vec[n]);
    const T_partials_return mu_dbl = value_of(mu_vec[n]);
    const T_partials_return lambda_dbl = value_of(lambda_vec[n]);
    const T_partials_return alpha_dbl = value_of(alpha_vec[n]);
    const T_partials_return temp = 1.0 + (y_dbl - mu_dbl) / lambda_dbl;

    const T_partials_return log_temp = log(temp);
    const T_partials_return rep_deriv
        = is_constant_all<T_y, T_loc, T_scale, T_shape>::value
              ? 0
              : alpha_dbl / (y_dbl - mu_dbl + lambda_dbl);

    const T_partials_return ccdf_log = -alpha_dbl * log_temp;

    P += ccdf_log;

    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_[n] -= rep_deriv;
    }
    if (!is_constant_all<T_loc>::value) {
      ops_partials.edge2_.partials_[n] += rep_deriv;
    }
    if (!is_constant_all<T_scale>::value) {
      ops_partials.edge3_.partials_[n]
          += rep_deriv * (y_dbl - mu_dbl) / lambda_dbl;
    }
    if (!is_constant_all<T_shape>::value) {
      ops_partials.edge4_.partials_[n] -= log_temp;
    }
  }

  return ops_partials.build(P);
}

}  // namespace math
}  // namespace stan
#endif
