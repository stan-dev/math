#ifndef STAN_MATH_PRIM_PROB_STD_NORMAL_LCDF_HPP
#define STAN_MATH_PRIM_PROB_STD_NORMAL_LCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/erf.hpp>
#include <stan/math/prim/fun/erfc.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1p.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>
#include <limits>

namespace stan {
namespace math {

template <
    typename T_y,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y>* = nullptr>
inline return_type_t<T_y> std_normal_lcdf(const T_y& y) {
  using T_partials_return = partials_return_t<T_y>;
  using std::exp;
  using std::fabs;
  using std::log;
  using std::pow;
  using T_y_ref = ref_type_t<T_y>;
  static const char* function = "std_normal_lcdf";
  T_y_ref y_ref = y;
  check_not_nan(function, "Random variable", y_ref);

  if (size_zero(y)) {
    return 0;
  }

  T_partials_return lcdf(0.0);
  auto ops_partials = make_partials_propagator(y_ref);

  scalar_seq_view<T_y_ref> y_vec(y_ref);
  size_t N = stan::math::size(y);

  for (size_t n = 0; n < N; n++) {
    const T_partials_return y_dbl = y_vec.val(n);
    const T_partials_return scaled_y = y_dbl * INV_SQRT_TWO;
    const T_partials_return x2 = square(scaled_y);

    // Rigorous numerical approximations are applied here to deal with values
    // of |scaled_y|>>0. This is needed to deal with rare base-rate
    // logistic regression problems where it is useful to use an alternative
    // link function instead.
    //
    // use erfc() instead of erf() in order to retain precision
    // since for x>0 erfc()->0
    if (scaled_y > 0.0) {
      // CDF(x) = 1/2 + 1/2erf(x) = 1 - 1/2erfc(x)
      lcdf += log1p(-0.5 * erfc(scaled_y));
      if (!is_not_nan(lcdf)) {
        lcdf = 0;
      }
    } else if (scaled_y > -20.0) {
      // CDF(x) = 1/2 - 1/2erf(-x) = 1/2erfc(-x)
      lcdf += log(erfc(-scaled_y)) + LOG_HALF;
    } else if (10.0 * log(fabs(scaled_y))
               < log(std::numeric_limits<T_partials_return>::max())) {
      // entering territory where erfc(-x)~0
      // need to use direct numerical approximation of lcdf instead
      // the following based on W. J. Cody, Math. Comp. 23(107):631-638 (1969)
      // CDF(x) = 1/2erfc(-x)
      const T_partials_return x4 = pow(scaled_y, 4);
      const T_partials_return x6 = pow(scaled_y, 6);
      const T_partials_return x8 = pow(scaled_y, 8);
      const T_partials_return x10 = pow(scaled_y, 10);
      const T_partials_return temp_p
          = 0.000658749161529837803157 + 0.0160837851487422766278 / x2
            + 0.125781726111229246204 / x4 + 0.360344899949804439429 / x6
            + 0.305326634961232344035 / x8 + 0.0163153871373020978498 / x10;
      const T_partials_return temp_q
          = -0.00233520497626869185443 - 0.0605183413124413191178 / x2
            - 0.527905102951428412248 / x4 - 1.87295284992346047209 / x6
            - 2.56852019228982242072 / x8 - 1.0 / x10;
      lcdf += LOG_HALF + log(INV_SQRT_PI + (temp_p / temp_q) / x2)
              - log(-scaled_y) - x2;
    } else {
      // scaled_y^10 term will overflow
      lcdf = stan::math::negative_infinity();
    }

    if (!is_constant_all<T_y>::value) {
      // compute partial derivatives
      // based on analytic form given by:
      // dln(CDF)/dx = exp(-x^2)/(sqrt(pi)*(1/2+erf(x)/2)
      T_partials_return dnlcdf = 0.0;
      T_partials_return t = 0.0;
      T_partials_return t2 = 0.0;
      T_partials_return t4 = 0.0;

      // calculate using piecewise function
      // (due to instability / inaccuracy in the various approximations)
      if (scaled_y > 2.9) {
        // approximation derived from Abramowitz and Stegun (1964) 7.1.26
        t = 1.0 / (1.0 + 0.3275911 * scaled_y);
        t2 = square(t);
        t4 = pow(t, 4);
        dnlcdf = INV_SQRT_PI
                 / (exp(x2) - 0.254829592 + 0.284496736 * t - 1.421413741 * t2
                    + 1.453152027 * t2 * t - 1.061405429 * t4);
      } else if (scaled_y > 2.5) {
        // in the trouble area where all of the standard numerical
        // approximations are unstable - bridge the gap using Taylor
        // expansions of the analytic function
        // use Taylor expansion centred around x=2.7
        t = scaled_y - 2.7;
        t2 = square(t);
        t4 = pow(t, 4);
        dnlcdf = 0.0003849882382 - 0.002079084702 * t + 0.005229340880 * t2
                 - 0.008029540137 * t2 * t + 0.008232190507 * t4
                 - 0.005692364250 * t4 * t + 0.002399496363 * pow(t, 6);
      } else if (scaled_y > 2.1) {
        // use Taylor expansion centred around x=2.3
        t = scaled_y - 2.3;
        t2 = square(t);
        t4 = pow(t, 4);
        dnlcdf = 0.002846135439 - 0.01310032351 * t + 0.02732189391 * t2
                 - 0.03326906904 * t2 * t + 0.02482478940 * t4
                 - 0.009883071924 * t4 * t - 0.0002771362254 * pow(t, 6);
      } else if (scaled_y > 1.5) {
        // use Taylor expansion centred around x=1.85
        t = scaled_y - 1.85;
        t2 = square(t);
        t4 = pow(t, 4);
        dnlcdf = 0.01849212058 - 0.06876280470 * t + 0.1099906382 * t2
                 - 0.09274533184 * t2 * t + 0.03543327418 * t4
                 + 0.005644855518 * t4 * t - 0.01111434424 * pow(t, 6);
      } else if (scaled_y > 0.8) {
        // use Taylor expansion centred around x=1.15
        t = scaled_y - 1.15;
        t2 = square(t);
        t4 = pow(t, 4);
        dnlcdf = 0.1585747034 - 0.3898677543 * t + 0.3515963775 * t2
                 - 0.09748053605 * t2 * t - 0.04347986191 * t4
                 + 0.02182506378 * t4 * t + 0.01074751427 * pow(t, 6);
      } else if (scaled_y > 0.1) {
        // use Taylor expansion centred around x=0.45
        t = scaled_y - 0.45;
        t2 = square(t);
        t4 = pow(t, 4);
        dnlcdf = 0.6245634904 - 0.9521866949 * t + 0.3986215682 * t2
                 + 0.04700850676 * t2 * t - 0.03478651979 * t4
                 - 0.01772675404 * t4 * t + 0.0006577254811 * pow(t, 6);
      } else if (10.0 * log(fabs(scaled_y))
                 < log(std::numeric_limits<T_partials_return>::max())) {
        // approximation derived from Abramowitz and Stegun (1964) 7.1.26
        // use fact that erf(x)=-erf(-x)
        // Abramowitz and Stegun define this for -inf<x<0 but seems to be
        // accurate for -inf<x<0.1
        t = 1.0 / (1.0 - 0.3275911 * scaled_y);
        t2 = square(t);
        t4 = pow(t, 4);
        dnlcdf = 2.0 * INV_SQRT_PI
                 / (0.254829592 * t - 0.284496736 * t2 + 1.421413741 * t2 * t
                    - 1.453152027 * t4 + 1.061405429 * t4 * t);
        // check if we need to add a correction term
        // (from cubic fit of residuals)
        if (scaled_y < -29.0) {
          dnlcdf += 0.0015065154280332 * x2 - 0.3993154819705530 * scaled_y
                    - 4.2919418242931700;
        } else if (scaled_y < -17.0) {
          dnlcdf += 0.0001263257217272 * x2 * scaled_y + 0.0123586859488623 * x2
                    - 0.0860505264736028 * scaled_y - 1.252783383752970;
        } else if (scaled_y < -7.0) {
          dnlcdf += 0.000471585349920831 * x2 * scaled_y
                    + 0.0296839305424034 * x2 + 0.207402143352332 * scaled_y
                    + 0.425316974683324;
        } else if (scaled_y < -3.9) {
          dnlcdf += -0.0006972280656443 * x2 * scaled_y
                    + 0.0068218494628567 * x2 + 0.0585761964460277 * scaled_y
                    + 0.1034397670201370;
        } else if (scaled_y < -2.1) {
          dnlcdf += -0.0018742199480885 * x2 * scaled_y
                    - 0.0097119598291202 * x2 - 0.0170137970924080 * scaled_y
                    - 0.0100428567412041;
        }
      } else {
        dnlcdf = stan::math::positive_infinity();
      }

      if (!is_constant_all<T_y>::value) {
        partials<0>(ops_partials)[n] += dnlcdf * INV_SQRT_TWO;
      }
    }
  }

  return ops_partials.build(lcdf);
}

}  // namespace math
}  // namespace stan
#endif
