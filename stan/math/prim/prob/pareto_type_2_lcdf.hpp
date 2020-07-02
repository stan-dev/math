#ifndef STAN_MATH_PRIM_PROB_PARETO_TYPE_2_LCDF_HPP
#define STAN_MATH_PRIM_PROB_PARETO_TYPE_2_LCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1m.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T_y, typename T_loc, typename T_scale, typename T_shape>
return_type_t<T_y, T_loc, T_scale, T_shape> pareto_type_2_lcdf(
    const T_y& y, const T_loc& mu, const T_scale& lambda,
    const T_shape& alpha) {
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale, T_shape>;
  using std::log;
  using std::pow;
  static const char* function = "pareto_type_2_lcdf";
  check_not_nan(function, "Random variable", y);
  check_nonnegative(function, "Random variable", y);
  check_positive_finite(function, "Scale parameter", lambda);
  check_positive_finite(function, "Shape parameter", alpha);
  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", lambda, "Shape parameter",
                         alpha);
  check_greater_or_equal(function, "Random variable", y, mu);

  if (size_zero(y, mu, lambda, alpha)) {
    return 0;
  }

  T_partials_return P(0.0);
  operands_and_partials<T_y, T_loc, T_scale, T_shape> ops_partials(
      y, mu, lambda, alpha);

  scalar_seq_view<T_y> y_vec(y);
  scalar_seq_view<T_loc> mu_vec(mu);
  scalar_seq_view<T_scale> lambda_vec(lambda);
  scalar_seq_view<T_shape> alpha_vec(alpha);
  size_t N = max_size(y, mu, lambda, alpha);

  for (size_t n = 0; n < N; n++) {
    const T_partials_return y_dbl = value_of(y_vec[n]);
    const T_partials_return mu_dbl = value_of(mu_vec[n]);
    const T_partials_return lambda_dbl = value_of(lambda_vec[n]);
    const T_partials_return alpha_dbl = value_of(alpha_vec[n]);
    const T_partials_return temp = 1.0 + (y_dbl - mu_dbl) / lambda_dbl;

    const T_partials_return p1_pow_alpha = pow(temp, alpha_dbl);
    const T_partials_return inv_p1_pow_alpha_minus_one
        = inv(p1_pow_alpha - 1.0);
    const T_partials_return grad_1_2
        = is_constant_all<T_y, T_loc, T_scale>::value
              ? 0
              : alpha_dbl * inv_p1_pow_alpha_minus_one
                    / (lambda_dbl - mu_dbl + y_dbl);

    const T_partials_return cdf_log = log1m(1.0 / p1_pow_alpha);

    P += cdf_log;

    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_[n] += grad_1_2;
    }
    if (!is_constant_all<T_loc>::value) {
      ops_partials.edge2_.partials_[n] -= grad_1_2;
    }
    if (!is_constant_all<T_scale>::value) {
      ops_partials.edge3_.partials_[n]
          += (mu_dbl - y_dbl) * grad_1_2 / lambda_dbl;
    }
    if (!is_constant_all<T_shape>::value) {
      ops_partials.edge4_.partials_[n]
          += log(temp) * inv_p1_pow_alpha_minus_one;
    }
  }

  return ops_partials.build(P);
}

}  // namespace math
}  // namespace stan
#endif
