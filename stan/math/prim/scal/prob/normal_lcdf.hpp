#ifndef STAN_MATH_PRIM_SCAL_PROB_NORMAL_LCDF_HPP
#define STAN_MATH_PRIM_SCAL_PROB_NORMAL_LCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/err/check_positive.hpp>
#include <stan/math/prim/scal/fun/size_zero.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <stan/math/prim/scal/fun/square.hpp>
#include <cmath>
#include <limits>

namespace stan {
namespace math {

template <typename T_y, typename T_loc, typename T_scale>
inline return_type_t<T_y, T_loc, T_scale> normal_lcdf(const T_y& y,
                                                      const T_loc& mu,
                                                      const T_scale& sigma) {
  static const char* function = "normal_lcdf";
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale>;
  using std::exp;
  using std::log;
  using std::log1p;
  using std::pow;

  T_partials_return cdf_log(0.0);
  if (size_zero(y, mu, sigma)) {
    return cdf_log;
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

  const T_partials_return SQRT_TWO_OVER_PI = std::sqrt(2.0 / pi());
  for (size_t n = 0; n < N; n++) {
    const T_partials_return y_dbl = value_of(y_vec[n]);
    const T_partials_return mu_dbl = value_of(mu_vec[n]);
    const T_partials_return sigma_dbl = value_of(sigma_vec[n]);

    const T_partials_return scaled_diff
        = (y_dbl - mu_dbl) / (sigma_dbl * SQRT_2);

    const T_partials_return sigma_sqrt2 = sigma_dbl * SQRT_2;
    const T_partials_return x2 = square(scaled_diff);

    // check whether scaled_diff^10 term will overflow
    const bool overflow_check
        = 10.0 * log(-scaled_diff)
          < log(std::numeric_limits<T_partials_return>::max());

    // use erfc() instead of erf() in order to retain precision
    // since for x>0 erfc()->0
    if (scaled_diff > 0.0) {
      // CDF(x) = 1/2 + 1/2erf(x) = 1 - 1/2erfc(x)
      cdf_log += log1p(-0.5 * erfc(scaled_diff));
      if (isnan(cdf_log)) {
        cdf_log = 0;
      }
    } else if (scaled_diff > -20.0) {
      // CDF(x) = 1/2 - 1/2erf(-x) = 1/2erfc(-x)
      cdf_log += log(erfc(-scaled_diff)) + LOG_HALF;
    } else {
      if (overflow_check) {
        // entering territory where erfc(-x)~0
        // need to use direct numerical approximation of cdf_log instead
        // the following based on W. J. Cody, Math. Comp. 23(107):631-638 (1969)
        // CDF(x) = 1/2erfc(-x)
        const T_partials_return x4 = pow(scaled_diff, 4);
        const T_partials_return x6 = pow(scaled_diff, 6);
        const T_partials_return x8 = pow(scaled_diff, 8);
        const T_partials_return x10 = pow(scaled_diff, 10);
        const T_partials_return temp_p
            = 0.000658749161529837803157 + 0.0160837851487422766278 / x2
              + 0.125781726111229246204 / x4 + 0.360344899949804439429 / x6
              + 0.305326634961232344035 / x8 + 0.0163153871373020978498 / x10;
        const T_partials_return temp_q
            = -0.00233520497626869185443 - 0.0605183413124413191178 / x2
              - 0.527905102951428412248 / x4 - 1.87295284992346047209 / x6
              - 2.56852019228982242072 / x8 - 1.0 / x10;
        cdf_log += LOG_HALF + log(1.0 / SQRT_PI + (temp_p / temp_q) / x2)
                   - log(-scaled_diff) - x2;
      } else {
        cdf_log = stan::math::negative_infinity();
      }
    }

    if (!is_constant_all<T_y, T_loc, T_scale>::value) {
      // compute partial derivatives
      // based on analytic form given by:
      // dln(CDF)/dx = exp(-x^2)/(sqrt(pi)*(1/2+erf(x)/2)
      T_partials_return dncdf_log_dbl = 0.0;
      T_partials_return t = 0.0, t2 = 0.0, t4 = 0.0;

      // calculate using piecewise funciton
      // (due to instability / inaccuracy in the various approximations)
      if (scaled_diff < 0.1) {
        // approximation derived from Abramowitx and Stegun (1964) 7.1.26
        // use fact that erf(x)=-erf(-x)
        // technically only true for -inf<x<0 but seems to be accurate
        // for -inf<x<0.1
        t = 1.0 / (1.0 - 0.3275911 * scaled_diff);
        t2 = square(t);
        t4 = pow(t, 4);
        if (overflow_check) {
          dncdf_log_dbl
              = 2.0
                / (SQRT_PI
                   * (0.254829592 * t - 0.284496736 * t2 + 1.421413741 * t2 * t
                      - 1.453152027 * t4 + 1.061405429 * t4 * t));
          if (scaled_diff < -4.0) {
            // need to add correction term (from cubic fit of residuals)
            dncdf_log_dbl += 0.00024073 * x2 * scaled_diff + 0.020656 * x2
                             + 0.10349 * scaled_diff + 0.97431;
          }
        } else {
          dncdf_log_dbl = stan::math::positive_infinity();
        }
      } else if (scaled_diff > 2.2) {
        // approximation derived from Abramowitx and Stegun (1964) 7.1.26
        t = 1.0 / (1.0 + 0.3275911 * scaled_diff);
        t2 = square(t);
        t4 = pow(t, 4);
        dncdf_log_dbl
            = 1.0
              / (SQRT_PI
                 * (exp(x2) - 0.254829592 + 0.284496736 * t - 1.421413741 * t2
                    + 1.453152027 * t2 * t - 1.061405429 * t4));
      } else {
        // in the trouble area where all of the standard numerical
        // approximations are unstable - bridge the gap using 3 Taylor
        // expansions of the analytic function
        if (scaled_diff < 0.8) {
          // use Taylor expansion centred around x=0.45
          t = scaled_diff - 0.45;
          t2 = square(t);
          t4 = pow(t, 4);
          dncdf_log_dbl = 0.6245634904 - 0.9521866949 * t + 0.3986215682 * t2
                          + 0.04700850676 * t2 * t - 0.03478651979 * t4
                          - 0.01772675404 * t4 * t
                          + 0.0006577254811 * pow(t, 6);
        } else if (scaled_diff < 1.5) {
          // use Taylor expansion centred around x=1.15
          t = scaled_diff - 1.15;
          t2 = square(t);
          t4 = pow(t, 4);
          dncdf_log_dbl = 0.1585747034 + 0.3898677543 * t + 0.3515963775 * t2
                          - 0.09748053605 * t2 * t - 0.04347986191 * t4
                          + 0.02182506378 * t4 * t + 0.01074751427 * pow(t, 6);
        } else {
          // use Taylor expansion centred around x=1.85
          t = scaled_diff - 1.85;
          t2 = square(t);
          t4 = pow(t, 4);
          dncdf_log_dbl = 0.01849212058 - 0.06876280470 * t + 0.1099906382 * t2
                          - 0.09274533184 * t2 * t + 0.03543327418 * t4
                          + 0.005644855518 * t4 * t - 0.01111434424 * pow(t, 6);
        }
      }

      if (!is_constant_all<T_y>::value) {
        ops_partials.edge1_.partials_[n] += dncdf_log_dbl / sigma_sqrt2;
      }
      if (!is_constant_all<T_loc>::value) {
        ops_partials.edge2_.partials_[n] -= dncdf_log_dbl / sigma_sqrt2;
      }
      if (!is_constant_all<T_scale>::value) {
        ops_partials.edge3_.partials_[n]
            -= dncdf_log_dbl * scaled_diff / sigma_dbl;
      }
    }
  }
  return ops_partials.build(cdf_log);
}

}  // namespace math
}  // namespace stan
#endif
