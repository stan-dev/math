#ifndef STAN_MATH_PRIM_SCAL_PROB_PARETO_TYPE_2_CDF_HPP
#define STAN_MATH_PRIM_SCAL_PROB_PARETO_TYPE_2_CDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_greater_or_equal.hpp>
#include <stan/math/prim/scal/err/check_nonnegative.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <stan/math/prim/scal/fun/size_zero.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T_y, typename T_loc, typename T_scale, typename T_shape>
return_type_t<T_y, T_loc, T_scale, T_shape> pareto_type_2_cdf(
    const T_y& y, const T_loc& mu, const T_scale& lambda,
    const T_shape& alpha) {
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale, T_shape>;

  if (size_zero(y, mu, lambda, alpha)) {
    return 1.0;
  }

  static const char* function = "pareto_type_2_cdf";

  using std::log;

  T_partials_return P(1.0);

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

  VectorBuilder<true, T_partials_return, T_y, T_loc, T_scale, T_shape>
      p1_pow_alpha(N);

  VectorBuilder<!is_constant_all<T_y, T_loc, T_scale>::value, T_partials_return,
                T_y, T_loc, T_scale, T_shape>
      grad_1_2(N);

  VectorBuilder<!is_constant_all<T_shape, T_y>::value, T_partials_return, T_y,
                T_loc, T_scale, T_shape>
      grad_3(N);

  for (size_t i = 0; i < N; i++) {
    const T_partials_return lambda_dbl = value_of(lambda_vec[i]);
    const T_partials_return alpha_dbl = value_of(alpha_vec[i]);
    const T_partials_return temp
        = 1 + (value_of(y_vec[i]) - value_of(mu_vec[i])) / lambda_dbl;
    p1_pow_alpha[i] = pow(temp, -alpha_dbl);

    if (!is_constant_all<T_y, T_loc, T_scale>::value) {
      grad_1_2[i] = p1_pow_alpha[i] / temp * alpha_dbl / lambda_dbl;
    }

    if (!is_constant_all<T_shape>::value) {
      grad_3[i] = log(temp) * p1_pow_alpha[i];
    }
  }

  for (size_t n = 0; n < N; n++) {
    const T_partials_return y_dbl = value_of(y_vec[n]);
    const T_partials_return mu_dbl = value_of(mu_vec[n]);
    const T_partials_return lambda_dbl = value_of(lambda_vec[n]);

    const T_partials_return Pn = 1.0 - p1_pow_alpha[n];

    P *= Pn;

    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_[n] += grad_1_2[n] / Pn;
    }
    if (!is_constant_all<T_loc>::value) {
      ops_partials.edge2_.partials_[n] -= grad_1_2[n] / Pn;
    }
    if (!is_constant_all<T_scale>::value) {
      ops_partials.edge3_.partials_[n]
          += (mu_dbl - y_dbl) * grad_1_2[n] / lambda_dbl / Pn;
    }
    if (!is_constant_all<T_shape>::value) {
      ops_partials.edge4_.partials_[n] += grad_3[n] / Pn;
    }
  }

  if (!is_constant_all<T_y>::value) {
    for (size_t n = 0; n < stan::length(y); ++n) {
      ops_partials.edge1_.partials_[n] *= P;
    }
  }
  if (!is_constant_all<T_loc>::value) {
    for (size_t n = 0; n < stan::length(mu); ++n) {
      ops_partials.edge2_.partials_[n] *= P;
    }
  }
  if (!is_constant_all<T_scale>::value) {
    for (size_t n = 0; n < stan::length(lambda); ++n) {
      ops_partials.edge3_.partials_[n] *= P;
    }
  }
  if (!is_constant_all<T_shape>::value) {
    for (size_t n = 0; n < stan::length(alpha); ++n) {
      ops_partials.edge4_.partials_[n] *= P;
    }
  }
  return ops_partials.build(P);
}

}  // namespace math
}  // namespace stan
#endif
