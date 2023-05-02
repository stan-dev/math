#ifndef STAN_MATH_PRIM_PROB_NORMAL_LCDF_HPP
#define STAN_MATH_PRIM_PROB_NORMAL_LCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/erf.hpp>
#include <stan/math/prim/fun/erfc.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1p.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>
#include <limits>

namespace stan {
namespace math {

template <typename T_y, typename T_loc, typename T_scale,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_loc, T_scale>* = nullptr>
inline return_type_t<T_y, T_loc, T_scale> normal_lcdf(const T_y& y,
                                                      const T_loc& mu,
                                                      const T_scale& sigma) {
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale>;
  using std::exp;
  using std::fabs;
  using std::log;
  using std::pow;
  using std::sqrt;
  using T_y_ref = ref_type_t<T_y>;
  using T_mu_ref = ref_type_t<T_loc>;
  using T_sigma_ref = ref_type_t<T_scale>;
  static const char* function = "normal_lcdf";
  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", sigma);
  T_y_ref y_ref = y;
  T_mu_ref mu_ref = mu;
  T_sigma_ref sigma_ref = sigma;
  check_not_nan(function, "Random variable", y_ref);
  check_finite(function, "Location parameter", mu_ref);
  check_positive(function, "Scale parameter", sigma_ref);

  if (size_zero(y, mu, sigma)) {
    return 0;
  }

  T_partials_return cdf_log(0.0);
  auto ops_partials = make_partials_propagator(y_ref, mu_ref, sigma_ref);

  scalar_seq_view<T_y_ref> y_vec(y_ref);
  scalar_seq_view<T_mu_ref> mu_vec(mu_ref);
  scalar_seq_view<T_sigma_ref> sigma_vec(sigma_ref);
  size_t N = max_size(y, mu, sigma);

  for (size_t n = 0; n < N; n++) {
    const T_partials_return y_dbl = y_vec.val(n);
    const T_partials_return mu_dbl = mu_vec.val(n);
    const T_partials_return sigma_dbl = sigma_vec.val(n);

    const T_partials_return scaled_diff
        = (y_dbl - mu_dbl) / (sigma_dbl * SQRT_TWO);

    const T_partials_return sigma_sqrt2 = sigma_dbl * SQRT_TWO;
    const T_partials_return x2 = square(scaled_diff);

    // Rigorous numerical approximations are applied here to deal with values
    // of |scaled_diff|>>0. This is needed to deal with rare base-rate
    // logistic regression problems where it is useful to use an alternative
    // link function instead.
    //
    // use erfc() instead of erf() in order to retain precision
    // since for x>0 erfc()->0
    if (scaled_diff > 0.0) {
      // CDF(x) = 1/2 + 1/2erf(x) = 1 - 1/2erfc(x)
      cdf_log += log1p(-0.5 * erfc(scaled_diff));
      if (!is_not_nan(cdf_log)) {
        cdf_log = 0;
      }
    } else if (scaled_diff > -20.0) {
      // CDF(x) = 1/2 - 1/2erf(-x) = 1/2erfc(-x)
      cdf_log += log(erfc(-scaled_diff)) + LOG_HALF;
    } else if (10.0 * log(fabs(scaled_diff))
               < log(std::numeric_limits<T_partials_return>::max())) {
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
      cdf_log += LOG_HALF + log(INV_SQRT_PI + (temp_p / temp_q) / x2)
                 - log(-scaled_diff) - x2;
    } else {
      // scaled_diff^10 term will overflow
      cdf_log = stan::math::negative_infinity();
    }

    if (!is_constant_all<T_y, T_loc, T_scale>::value) {
      // compute partial derivatives
      // based on analytic form given by:
      // dln(CDF)/dx = exp(-x^2)/(sqrt(pi)*(1/2+erf(x)/2)
      T_partials_return dncdf_log = 0.0;
      T_partials_return t = 0.0;
      T_partials_return t2 = 0.0;
      T_partials_return t4 = 0.0;

      // calculate using piecewise function
      // (due to instability / inaccuracy in the various approximations)
      if (scaled_diff > 2.9) {
        // approximation derived from Abramowitz and Stegun (1964) 7.1.26
        t = 1.0 / (1.0 + 0.3275911 * scaled_diff);
        t2 = square(t);
        t4 = pow(t, 4);
        dncdf_log
            = 1.0
              / (SQRT_PI
                 * (exp(x2) - 0.254829592 + 0.284496736 * t - 1.421413741 * t2
                    + 1.453152027 * t2 * t - 1.061405429 * t4));
      } else if (scaled_diff > 2.5) {
        // in the trouble area where all of the standard numerical
        // approximations are unstable - bridge the gap using Taylor
        // expansions of the analytic function
        // use Taylor expansion centred around x=2.7
        t = scaled_diff - 2.7;
        t2 = square(t);
        t4 = pow(t, 4);
        dncdf_log = 0.0003849882382 - 0.002079084702 * t + 0.005229340880 * t2
                    - 0.008029540137 * t2 * t + 0.008232190507 * t4
                    - 0.005692364250 * t4 * t + 0.002399496363 * pow(t, 6);
      } else if (scaled_diff > 2.1) {
        // use Taylor expansion centred around x=2.3
        t = scaled_diff - 2.3;
        t2 = square(t);
        t4 = pow(t, 4);
        dncdf_log = 0.002846135439 - 0.01310032351 * t + 0.02732189391 * t2
                    - 0.03326906904 * t2 * t + 0.02482478940 * t4
                    - 0.009883071924 * t4 * t - 0.0002771362254 * pow(t, 6);
      } else if (scaled_diff > 1.5) {
        // use Taylor expansion centred around x=1.85
        t = scaled_diff - 1.85;
        t2 = square(t);
        t4 = pow(t, 4);
        dncdf_log = 0.01849212058 - 0.06876280470 * t + 0.1099906382 * t2
                    - 0.09274533184 * t2 * t + 0.03543327418 * t4
                    + 0.005644855518 * t4 * t - 0.01111434424 * pow(t, 6);
      } else if (scaled_diff > 0.8) {
        // use Taylor expansion centred around x=1.15
        t = scaled_diff - 1.15;
        t2 = square(t);
        t4 = pow(t, 4);
        dncdf_log = 0.1585747034 - 0.3898677543 * t + 0.3515963775 * t2
                    - 0.09748053605 * t2 * t - 0.04347986191 * t4
                    + 0.02182506378 * t4 * t + 0.01074751427 * pow(t, 6);
      } else if (scaled_diff > 0.1) {
        // use Taylor expansion centred around x=0.45
        t = scaled_diff - 0.45;
        t2 = square(t);
        t4 = pow(t, 4);
        dncdf_log = 0.6245634904 - 0.9521866949 * t + 0.3986215682 * t2
                    + 0.04700850676 * t2 * t - 0.03478651979 * t4
                    - 0.01772675404 * t4 * t + 0.0006577254811 * pow(t, 6);
      } else if (10.0 * log(fabs(scaled_diff))
                 < log(std::numeric_limits<T_partials_return>::max())) {
        // approximation derived from Abramowitz and Stegun (1964) 7.1.26
        // use fact that erf(x)=-erf(-x)
        // Abramowitz and Stegun define this for -inf<x<0 but seems to be
        // accurate for -inf<x<0.1
        t = 1.0 / (1.0 - 0.3275911 * scaled_diff);
        t2 = square(t);
        t4 = pow(t, 4);
        dncdf_log
            = 2.0
              / (SQRT_PI
                 * (0.254829592 * t - 0.284496736 * t2 + 1.421413741 * t2 * t
                    - 1.453152027 * t4 + 1.061405429 * t4 * t));
        // check if we need to add a correction term
        // (from cubic fit of residuals)
        if (scaled_diff < -29.0) {
          dncdf_log += 0.0015065154280332 * x2
                       - 0.3993154819705530 * scaled_diff - 4.2919418242931700;
        } else if (scaled_diff < -17.0) {
          dncdf_log += 0.0001263257217272 * x2 * scaled_diff
                       + 0.0123586859488623 * x2
                       - 0.0860505264736028 * scaled_diff - 1.252783383752970;
        } else if (scaled_diff < -7.0) {
          dncdf_log += 0.000471585349920831 * x2 * scaled_diff
                       + 0.0296839305424034 * x2
                       + 0.207402143352332 * scaled_diff + 0.425316974683324;
        } else if (scaled_diff < -3.9) {
          dncdf_log += -0.0006972280656443 * x2 * scaled_diff
                       + 0.0068218494628567 * x2
                       + 0.0585761964460277 * scaled_diff + 0.1034397670201370;
        } else if (scaled_diff < -2.1) {
          dncdf_log += -0.0018742199480885 * x2 * scaled_diff
                       - 0.0097119598291202 * x2
                       - 0.0170137970924080 * scaled_diff - 0.0100428567412041;
        }
      } else {
        dncdf_log = stan::math::positive_infinity();
      }

      if (!is_constant_all<T_y>::value) {
        partials<0>(ops_partials)[n] += dncdf_log / sigma_sqrt2;
      }
      if (!is_constant_all<T_loc>::value) {
        partials<1>(ops_partials)[n] -= dncdf_log / sigma_sqrt2;
      }
      if (!is_constant_all<T_scale>::value) {
        partials<2>(ops_partials)[n] -= dncdf_log * scaled_diff / sigma_dbl;
      }
    }
  }
  return ops_partials.build(cdf_log);
}

}  // namespace math
}  // namespace stan
#endif
