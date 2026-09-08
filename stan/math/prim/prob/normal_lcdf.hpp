#ifndef STAN_MATH_PRIM_PROB_NORMAL_LCDF_HPP
#define STAN_MATH_PRIM_PROB_NORMAL_LCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/erf.hpp>
#include <stan/math/prim/fun/erfc.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/fabs.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1p.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/pow.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/fun/sqrt.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>
#include <limits>

namespace stan {
namespace math {
namespace internal {
constexpr char normal_lcdf_func[] = "normal_lcdf";
}  // namespace internal

/** \ingroup prob_dists
 * @brief Calculates the log of the cdf of the normal distribution
 *
 * Tail branching follows three published results, and matches what other
 * libraries do:
 *  - Abramowitz & Stegun (1964) 7.1.26 Page 299 gives
 *    `erf(x) = 1 - P(t) exp(-x^2)`, `t = 1/(1 + p x)`, with `exp(-x^2)` as a
 *    numerator factor:
 *    https://archive.org/details/handbookofmathem1964abra/page/298/mode/2up
 *  - Cody (1969), Math. Comp. 23:631-637, is the rational approximation used
 *    for the tail value: https://doi.org/10.1090/S0025-5718-1969-0247736-4
 *  - Digital Library of Mathematical Functions (DLMF) 7.12.1
 *    is the erfc asymptotic behind the far-tail Mills ratio:
 *    https://dlmf.nist.gov/7.12.E1
 *
 * R's `pnorm` only ever forms the Gaussian factor in the numerator
 * (the `do_del` macro), and switches to the Cody tail form at `y > M_SQRT_32`,
 * i.e. `|x| > sqrt(32) ~= 5.657`. Since `scaled_diff = x / sqrt(2)`, that same
 * crossover is `|scaled_diff| > sqrt(32)/sqrt(2) = 4` exactly, which is the
 * cutoff used here for the cdf value. Measured by evaluating the
 * `temp_p`/`temp_q` expression below against `log(erfc(-scaled_diff)/2)` at 60
 * significant digits, our Cody set holds to 8.6e-17 relative down to
 * scaled_diff = -4 and degrades past about -3.52, so 4 sits just inside its
 * range. This is the same as R's impl.
 * https://github.com/wch/r-source/blob/trunk/src/nmath/pnorm.c
 * SciPy's `log_ndtr` uses the identical `log1p(-erfc(t)/2)` upper branch with
 * the same `t = x/sqrt(2)`, and needs no rational approximation or Taylor
 * patches at all below `x = -1`, because `erfcx` never forms `exp(+t^2)`:
 * https://github.com/scipy/xsf/blob/main/include/xsf/stats.h
 *
 * The interior cutoffs for the gradient are not from the literature.
 * They were derived in stan-dev/math#1411
 * (Phil Clemson, Univ. of Liverpool, Nov 2019), fixing
 * stan-dev/math#1284, which describes them as "original Taylor expansions
 * that have been derived to bridge the gap where the numerical
 * approximations are unstable". The author's account of how they were
 * placed, from the review thread, was: "After playing around with the
 * autodiff tester I found some regions where it was failing some of the
 * tests, so I added some new approximations (Taylor expansions and fits of
 * the residuals) to improve the accuracy."
 * https://github.com/stan-dev/math/pull/1411
 * https://github.com/stan-dev/math/issues/1284
 *
 * Each row of the chart below is one Taylor expansion around `scaled_diff`.
 * `interval` is the range of `scaled_diff` it covers.
 * `centre` is the point the series is expanded about.
 *
 * | interval | centre | half-width |
 * |---|---|---|
 * | (2.5, 2.9] | 2.7 | 0.200 |
 * | (2.1, 2.5] | 2.3 | 0.200 |
 * | (1.5, 2.1] | 1.85 | 0.300 |
 * | (0.8, 1.5] | 1.15 | 0.350 |
 * | (0.1, 0.8] | 0.45 | 0.350 |
 *
 * The cut points are not arbitrary. Each interval of `scaled_diff` is centred
 * on its own Taylor point, which minimises the largest `|t|` the series has to
 * cover.
 *
 * The second table below shows why each interval ends where it does: the
 * series is accurate across its own range and falls apart just outside it,
 * so the intervals cannot be widened to use fewer of them.
 *
 * Worst in-range comes from evaluating each Taylor polynomial as written
 * against `(2/sqrt(pi)) * exp(-scaled_diff^2) / erfc(-scaled_diff)`, the exact
 * derivative, at 60 significant digits.
 *
 * Columns below:
 * `worst in-range`: The largest relative error of that branch's series over
 *   its own interval;
 * `0.2 below lo`: The relative error the same series would give at `lo - 0.2`;
 * `0.2 above hi`: The relative error the same series would give at `hi + 0.2`,
 *   i.e. what widening the interval either way would cost.
 *
 * | interval | worst in-range | 0.2 below lo | 0.2 above hi |
 * |---|---|---|---|
 * | (2.5, 2.9] | 3.27e-05 | 2.32e-05 | 1.61e-02 |
 * | (2.1, 2.5] | 2.80e-05 | 3.25e-04 | 9.15e-03 |
 * | (1.5, 2.1] | 2.30e-05 | 8.11e-05 | 3.77e-03 |
 * | (0.8, 1.5] | 6.09e-05 | 9.29e-05 | 2.93e-03 |
 * | (0.1, 0.8] | 7.61e-06 | 3.62e-05 | 2.82e-04 |
 *
 * The worst in-range error is uniformly 1e-5 to 6e-5, just inside the 1e-4
 * relative gradient tolerance `expect_ad` applies by default
 * (`gradient_grad_` in test/unit/math/ad_tolerances.hpp; the 1e-8 there is
 * `gradient_val_`, which bounds the value rather than the gradient), while
 * widening an interval by 0.2 costs one to three orders of magnitude. That
 * uniformity, not any published result, is the placement criterion.
 *
 * The `worst in-range` column is enforced: the branch_accuracy test in
 * mix/prob/normal_cdf_log_test.cpp asserts every branch against
 * 60-digit references at that error plus a small margin.
 *
 * The negative-tail residual corrections at `scaled_diff` of -2.1, -3.9, -7 and
 * -17 are cubic fits of residuals from the same PR and likewise have no
 * external source. The -29 cutoff below is justified by DLMF 7.12.1
 *
 * @tparam func name reported by the error checks. Reflected distributions
 *   such as `normal_lccdf` delegate here and pass their own name so that
 *   exceptions name the function the user actually called.
 * @tparam T_y A vector or scalar type for the random variable.
 * @tparam T_loc A vector or scalar type for the location parameter.
 * @tparam T_scale A vector or scalar type for the scale parameter.
 * @param y (Sequence of) scalar(s).
 * @param mu (Sequence of) scalar(s).
 * @param sigma (Sequence of) scalar(s).
 * @return The log of the normal cdf evaluated at the specified arguments. If
 *   given containers, the log of the product of the cdfs.
 */
template <const char* func = internal::normal_lcdf_func, typename T_y,
          typename T_loc, typename T_scale,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_loc, T_scale>* = nullptr>
inline return_type_t<T_y, T_loc, T_scale> normal_lcdf(const T_y& y,
                                                      const T_loc& mu,
                                                      const T_scale& sigma) {
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale>;
  using T_y_ref = ref_type_t<T_y>;
  using T_mu_ref = ref_type_t<T_loc>;
  using T_sigma_ref = ref_type_t<T_scale>;
  static constexpr const char* function = func;
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
    } else if (scaled_diff > -4.0) {
      // CDF(x) = 1/2 - 1/2erf(-x) = 1/2erfc(-x); -4 is R pnorm's M_SQRT_32
      // Since we scale by sqrt(2), we use sqrt(32)/sqrt(2) = 4
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

    if constexpr (is_any_autodiff_v<T_y, T_loc, T_scale>) {
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
        // A&S 7.1.26 puts exp(-x2) in the numerator; keep it there so it
        // underflows to zero instead of overflowing inside a denominator
        const T_partials_return exp_m_x2 = exp(-x2);
        dncdf_log
            = exp_m_x2
              / (SQRT_PI
                 * (1.0
                    - exp_m_x2
                          * (0.254829592 - 0.284496736 * t + 1.421413741 * t2
                             - 1.453152027 * t2 * t + 1.061405429 * t4)));
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
      } else if (scaled_diff < -29.0) {
        // asymptotic Mills ratio, DLMF 7.12.1: dncdf_log grows linearly as
        // -2*scaled_diff, so no quadratic residual fit can track it
        const T_partials_return inv_x2 = 1.0 / x2;
        dncdf_log
            = -2.0 * scaled_diff
              / (1.0 + inv_x2 * (-0.5 + inv_x2 * (0.75 + inv_x2 * -1.875)));
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
        if (scaled_diff < -17.0) {
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
      const T_partials_return sigma_sqrt2 = sigma_dbl * SQRT_TWO;
      if constexpr (is_autodiff_v<T_y>) {
        partials<0>(ops_partials)[n] += dncdf_log / sigma_sqrt2;
      }
      if constexpr (is_autodiff_v<T_loc>) {
        partials<1>(ops_partials)[n] -= dncdf_log / sigma_sqrt2;
      }
      if constexpr (is_autodiff_v<T_scale>) {
        partials<2>(ops_partials)[n] -= dncdf_log * scaled_diff / sigma_dbl;
      }
    }
  }
  return ops_partials.build(cdf_log);
}

}  // namespace math
}  // namespace stan
#endif
