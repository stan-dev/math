#ifndef STAN_MATH_PRIM_PROB_NORMAL_LCCDF_HPP
#define STAN_MATH_PRIM_PROB_NORMAL_LCCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/erf.hpp>
#include <stan/math/prim/fun/erfc.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T_y, typename T_loc, typename T_scale>
inline return_type_t<T_y, T_loc, T_scale> normal_lccdf(const T_y& y,
                                                       const T_loc& mu,
                                                       const T_scale& sigma) {
  static const char* function = "normal_lccdf";
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale>;

  using std::exp;
  using std::log;

  T_partials_return ccdf_log(0.0);
  if (size_zero(y, mu, sigma)) {
    return ccdf_log;
  }

  check_not_nan(function, "Random variable", y);
  check_finite(function, "Location parameter", mu);
  check_not_nan(function, "Scale parameter", sigma);
  check_positive(function, "Scale parameter", sigma);
  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", sigma);

  operands_and_partials<T_y, T_loc, T_scale> ops_partials(y, mu, sigma);

  scalar_seq_view<T_y> y_vec(y);
  scalar_seq_view<T_loc> mu_vec(mu);
  scalar_seq_view<T_scale> sigma_vec(sigma);
  size_t N = max_size(y, mu, sigma);

  for (size_t n = 0; n < N; n++) {
    const T_partials_return y_dbl = value_of(y_vec[n]);
    const T_partials_return mu_dbl = value_of(mu_vec[n]);
    const T_partials_return sigma_dbl = value_of(sigma_vec[n]);

    const T_partials_return scaled_diff
        = (y_dbl - mu_dbl) / (sigma_dbl * SQRT_TWO);

    T_partials_return one_m_erf;
    if (scaled_diff < -37.5 * INV_SQRT_TWO) {
      one_m_erf = 2.0;
    } else if (scaled_diff < -5.0 * INV_SQRT_TWO) {
      one_m_erf = 2.0 - erfc(-scaled_diff);
    } else if (scaled_diff > 8.25 * INV_SQRT_TWO) {
      one_m_erf = 0.0;
    } else {
      one_m_erf = 1.0 - erf(scaled_diff);
    }

    ccdf_log += LOG_HALF + log(one_m_erf);

    if (!is_constant_all<T_y, T_loc, T_scale>::value) {
      const T_partials_return rep_deriv_div_sigma
          = scaled_diff > 8.25 * INV_SQRT_TWO
                ? INFTY
                : SQRT_TWO_OVER_SQRT_PI * exp(-scaled_diff * scaled_diff)
                      / one_m_erf / sigma_dbl;
      if (!is_constant_all<T_y>::value) {
        ops_partials.edge1_.partials_[n] -= rep_deriv_div_sigma;
      }
      if (!is_constant_all<T_loc>::value) {
        ops_partials.edge2_.partials_[n] += rep_deriv_div_sigma;
      }
      if (!is_constant_all<T_scale>::value) {
        ops_partials.edge3_.partials_[n]
            += rep_deriv_div_sigma * scaled_diff * SQRT_TWO;
      }
    }
  }
  return ops_partials.build(ccdf_log);
}

}  // namespace math
}  // namespace stan
#endif
