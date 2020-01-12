#ifndef STAN_MATH_PRIM_PROB_LOGNORMAL_LCDF_HPP
#define STAN_MATH_PRIM_PROB_LOGNORMAL_LCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/erfc.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T_y, typename T_loc, typename T_scale>
return_type_t<T_y, T_loc, T_scale> lognormal_lcdf(const T_y& y, const T_loc& mu,
                                                  const T_scale& sigma) {
  static const char* function = "lognormal_lcdf";
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale>;

  T_partials_return cdf_log = 0.0;

  using std::exp;
  using std::log;

  if (size_zero(y, mu, sigma)) {
    return cdf_log;
  }

  check_not_nan(function, "Random variable", y);
  check_nonnegative(function, "Random variable", y);
  check_finite(function, "Location parameter", mu);
  check_positive_finite(function, "Scale parameter", sigma);

  operands_and_partials<T_y, T_loc, T_scale> ops_partials(y, mu, sigma);

  scalar_seq_view<T_y> y_vec(y);
  scalar_seq_view<T_loc> mu_vec(mu);
  scalar_seq_view<T_scale> sigma_vec(sigma);
  size_t N = max_size(y, mu, sigma);

  for (size_t i = 0; i < size(y); i++) {
    if (value_of(y_vec[i]) == 0.0) {
      return ops_partials.build(negative_infinity());
    }
  }

  for (size_t n = 0; n < N; n++) {
    const T_partials_return y_dbl = value_of(y_vec[n]);
    const T_partials_return mu_dbl = value_of(mu_vec[n]);
    const T_partials_return sigma_dbl = value_of(sigma_vec[n]);
    const T_partials_return scaled_diff
        = (log(y_dbl) - mu_dbl) / (sigma_dbl * SQRT_TWO);
    const T_partials_return rep_deriv
        = SQRT_TWO_OVER_SQRT_PI * exp(-scaled_diff * scaled_diff) / sigma_dbl;

    const T_partials_return erfc_calc = erfc(-scaled_diff);
    cdf_log += LOG_HALF + log(erfc_calc);

    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_[n] += rep_deriv / erfc_calc / y_dbl;
    }
    if (!is_constant_all<T_loc>::value) {
      ops_partials.edge2_.partials_[n] -= rep_deriv / erfc_calc;
    }
    if (!is_constant_all<T_scale>::value) {
      ops_partials.edge3_.partials_[n]
          -= rep_deriv * scaled_diff * SQRT_TWO / erfc_calc;
    }
  }
  return ops_partials.build(cdf_log);
}

}  // namespace math
}  // namespace stan
#endif
