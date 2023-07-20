#ifndef STAN_MATH_OPENCL_PRIM_NORMAL_LCDF_HPP
#define STAN_MATH_OPENCL_PRIM_NORMAL_LCDF_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/elt_divide.hpp>
#include <stan/math/prim/fun/elt_multiply.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>

namespace stan {
namespace math {
namespace internal {
const char opencl_normal_lcdf_impl[] = STRINGIFY(
    double x2 = normal_lcdf_scaled_diff * normal_lcdf_scaled_diff;
    double normal_lcdf_n = 0;
    // Rigorous numerical approximations are applied here to deal with values
    // of |normal_lcdf_scaled_diff|>>0. This is needed to deal with rare
    // base-rate logistic regression problems where it is useful to use an
    // alternative link function instead.
    //
    // use erfc() instead of erf() in order to retain precision
    // since for x>0 erfc()->0
    if (normal_lcdf_scaled_diff > 0.0) {
      // CDF(x) = 1/2 + 1/2erf(x) = 1 - 1/2erfc(x)
      normal_lcdf_n = log1p(-0.5 * erfc(normal_lcdf_scaled_diff));
      if (isnan(normal_lcdf_n)) {
        normal_lcdf_n = 0;
      }
    } else if (normal_lcdf_scaled_diff > -20.0) {
      // CDF(x) = 1/2 - 1/2erf(-x) = 1/2erfc(-x)
      normal_lcdf_n = log(erfc(-normal_lcdf_scaled_diff)) - M_LN2;
    } else if (10.0 * log(fabs(normal_lcdf_scaled_diff)) < log(DBL_MAX)) {
      // entering territory where erfc(-x)~0
      // need to use direct numerical approximation of normal_lcdf_n instead
      // the following based on W. J. Cody, Math. Comp. 23(107):631-638 (1969)
      // CDF(x) = 1/2erfc(-x)
      double x4 = pow(normal_lcdf_scaled_diff, 4);
      double x6 = pow(normal_lcdf_scaled_diff, 6);
      double x8 = pow(normal_lcdf_scaled_diff, 8);
      double x10 = pow(normal_lcdf_scaled_diff, 10);
      double temp_p
          = 0.000658749161529837803157 + 0.0160837851487422766278 / x2
            + 0.125781726111229246204 / x4 + 0.360344899949804439429 / x6
            + 0.305326634961232344035 / x8 + 0.0163153871373020978498 / x10;
      double temp_q = -0.00233520497626869185443 - 0.0605183413124413191178 / x2
                      - 0.527905102951428412248 / x4
                      - 1.87295284992346047209 / x6
                      - 2.56852019228982242072 / x8 - 1.0 / x10;
      normal_lcdf_n = -M_LN2 + log(0.5 * M_2_SQRTPI + (temp_p / temp_q) / x2)
                      - log(-normal_lcdf_scaled_diff) - x2;
    } else {
      // normal_lcdf_scaled_diff^10 term will overflow
      normal_lcdf_n = -INFINITY;
    });

const char opencl_normal_lcdf_ldncdf_impl[] = STRINGIFY(
    double normal_ldncdf = 0.0; double t = 0.0; double t2 = 0.0;
    double t4 = 0.0;

    // calculate using piecewise function
    // (due to instability / inaccuracy in the various approximations)
    if (normal_lcdf_deriv_scaled_diff > 2.9) {
      // approximation derived from Abramowitz and Stegun (1964) 7.1.26
      t = 1.0 / (1.0 + 0.3275911 * normal_lcdf_deriv_scaled_diff);
      t2 = t * t;
      t4 = pow(t, 4);
      normal_ldncdf
          = 0.5 * M_2_SQRTPI
            / (exp(x2) - 0.254829592 + 0.284496736 * t - 1.421413741 * t2
               + 1.453152027 * t2 * t - 1.061405429 * t4);
    } else if (normal_lcdf_deriv_scaled_diff > 2.5) {
      // in the trouble area where all of the standard numerical
      // approximations are unstable - bridge the gap using Taylor
      // expansions of the analytic function
      // use Taylor expansion centred around x=2.7
      t = normal_lcdf_deriv_scaled_diff - 2.7;
      t2 = t * t;
      t4 = pow(t, 4);
      normal_ldncdf = 0.0003849882382 - 0.002079084702 * t + 0.005229340880 * t2
                      - 0.008029540137 * t2 * t + 0.008232190507 * t4
                      - 0.005692364250 * t4 * t + 0.002399496363 * pow(t, 6);
    } else if (normal_lcdf_deriv_scaled_diff > 2.1) {
      // use Taylor expansion centred around x=2.3
      t = normal_lcdf_deriv_scaled_diff - 2.3;
      t2 = t * t;
      t4 = pow(t, 4);
      normal_ldncdf = 0.002846135439 - 0.01310032351 * t + 0.02732189391 * t2
                      - 0.03326906904 * t2 * t + 0.02482478940 * t4
                      - 0.009883071924 * t4 * t - 0.0002771362254 * pow(t, 6);
    } else if (normal_lcdf_deriv_scaled_diff > 1.5) {
      // use Taylor expansion centred around x=1.85
      t = normal_lcdf_deriv_scaled_diff - 1.85;
      t2 = t * t;
      t4 = pow(t, 4);
      normal_ldncdf = 0.01849212058 - 0.06876280470 * t + 0.1099906382 * t2
                      - 0.09274533184 * t2 * t + 0.03543327418 * t4
                      + 0.005644855518 * t4 * t - 0.01111434424 * pow(t, 6);
    } else if (normal_lcdf_deriv_scaled_diff > 0.8) {
      // use Taylor expansion centred around x=1.15
      t = normal_lcdf_deriv_scaled_diff - 1.15;
      t2 = t * t;
      t4 = pow(t, 4);
      normal_ldncdf = 0.1585747034 - 0.3898677543 * t + 0.3515963775 * t2
                      - 0.09748053605 * t2 * t - 0.04347986191 * t4
                      + 0.02182506378 * t4 * t + 0.01074751427 * pow(t, 6);
    } else if (normal_lcdf_deriv_scaled_diff > 0.1) {
      // use Taylor expansion centred around x=0.45
      t = normal_lcdf_deriv_scaled_diff - 0.45;
      t2 = t * t;
      t4 = pow(t, 4);
      normal_ldncdf = 0.6245634904 - 0.9521866949 * t + 0.3986215682 * t2
                      + 0.04700850676 * t2 * t - 0.03478651979 * t4
                      - 0.01772675404 * t4 * t + 0.0006577254811 * pow(t, 6);
    } else if (10.0 * log(fabs(normal_lcdf_deriv_scaled_diff)) < log(DBL_MAX)) {
      // approximation derived from Abramowitz and Stegun (1964) 7.1.26
      // use fact that erf(x)=-erf(-x)
      // Abramowitz and Stegun define this for -inf<x<0 but seems to be
      // accurate for -inf<x<0.1
      t = 1.0 / (1.0 - 0.3275911 * normal_lcdf_deriv_scaled_diff);
      t2 = t * t;
      t4 = pow(t, 4);
      normal_ldncdf
          = M_2_SQRTPI
            / (0.254829592 * t - 0.284496736 * t2 + 1.421413741 * t2 * t
               - 1.453152027 * t4 + 1.061405429 * t4 * t);
      // check if we need to add a correction term
      // (from cubic fit of residuals)
      if (normal_lcdf_deriv_scaled_diff < -29.0) {
        normal_ldncdf += 0.0015065154280332 * x2
                         - 0.3993154819705530 * normal_lcdf_deriv_scaled_diff
                         - 4.2919418242931700;
      } else if (normal_lcdf_deriv_scaled_diff < -17.0) {
        normal_ldncdf += 0.0001263257217272 * x2 * normal_lcdf_deriv_scaled_diff
                         + 0.0123586859488623 * x2
                         - 0.0860505264736028 * normal_lcdf_deriv_scaled_diff
                         - 1.252783383752970;
      } else if (normal_lcdf_deriv_scaled_diff < -7.0) {
        normal_ldncdf
            += 0.000471585349920831 * x2 * normal_lcdf_deriv_scaled_diff
               + 0.0296839305424034 * x2
               + 0.207402143352332 * normal_lcdf_deriv_scaled_diff
               + 0.425316974683324;
      } else if (normal_lcdf_deriv_scaled_diff < -3.9) {
        normal_ldncdf
            += -0.0006972280656443 * x2 * normal_lcdf_deriv_scaled_diff
               + 0.0068218494628567 * x2
               + 0.0585761964460277 * normal_lcdf_deriv_scaled_diff
               + 0.1034397670201370;
      } else if (normal_lcdf_deriv_scaled_diff < -2.1) {
        normal_ldncdf
            += -0.0018742199480885 * x2 * normal_lcdf_deriv_scaled_diff
               - 0.0097119598291202 * x2
               - 0.0170137970924080 * normal_lcdf_deriv_scaled_diff
               - 0.0100428567412041;
      }
    } else { normal_ldncdf = INFINITY; });
}  // namespace internal

/** \ingroup opencl
 * Returns the normal log complementary cumulative distribution function
 * for the given location, and scale. If given containers of matching sizes
 * returns the log sum of probabilities.
 *
 * @tparam T_y_cl type of scalar outcome
 * @tparam T_loc_cl type of location
 * @tparam T_scale_cl type of scale
 * @param y (Sequence of) scalar(s).
 * @param mu (Sequence of) location(s).
 * @param sigma (Sequence of) scale(s).
 * @return The log of the product of densities.
 */
template <
    typename T_y_cl, typename T_loc_cl, typename T_scale_cl,
    require_all_prim_or_rev_kernel_expression_t<T_y_cl, T_loc_cl,
                                                T_scale_cl>* = nullptr,
    require_any_not_stan_scalar_t<T_y_cl, T_loc_cl, T_scale_cl>* = nullptr>
return_type_t<T_y_cl, T_loc_cl, T_scale_cl> normal_lcdf(
    const T_y_cl& y, const T_loc_cl& mu, const T_scale_cl& sigma) {
  static constexpr const char* function = "normal_lcdf(OpenCL)";
  using std::isfinite;
  using std::isnan;

  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", sigma);
  const size_t N = max_size(y, mu, sigma);
  if (N == 0) {
    return 0.0;
  }

  const auto& y_col = as_column_vector_or_scalar(y);
  const auto& mu_col = as_column_vector_or_scalar(mu);
  const auto& sigma_col = as_column_vector_or_scalar(sigma);

  const auto& y_val = value_of(y_col);
  const auto& mu_val = value_of(mu_col);
  const auto& sigma_val = value_of(sigma_col);

  auto check_y_not_nan
      = check_cl(function, "Random variable", y_val, "not NaN");
  auto y_not_nan_expr = !isnan(y_val);
  auto check_mu_finite
      = check_cl(function, "Location parameter", mu_val, "finite");
  auto mu_finite_expr = isfinite(mu_val);
  auto check_sigma_positive
      = check_cl(function, "Scale parameter", sigma_val, "positive");
  auto sigma_positive_expr = 0 < sigma_val;

  auto scaled_diff = elt_divide(y_val - mu_val, sigma_val * SQRT_TWO);

  auto sigma_sqrt2 = sigma_val * SQRT_TWO;

  auto lcdf_n = opencl_code<internal::opencl_normal_lcdf_impl>(
                    std::make_tuple("normal_lcdf_scaled_diff"), scaled_diff)
                    .template output<double>("normal_lcdf_n");
  auto lcdf_expr = colwise_sum(lcdf_n);

  auto ldncdf
      = opencl_code<internal::opencl_normal_lcdf_ldncdf_impl>(
            std::make_tuple("normal_lcdf_deriv_scaled_diff"), scaled_diff)
            .template output<double>("normal_ldncdf");
  auto y_deriv = elt_divide(ldncdf, sigma_sqrt2);
  auto mu_deriv = -y_deriv;
  auto sigma_deriv = -elt_divide(elt_multiply(ldncdf, scaled_diff), sigma_val);

  matrix_cl<double> lcdf_cl;
  matrix_cl<double> y_deriv_cl;
  matrix_cl<double> mu_deriv_cl;
  matrix_cl<double> sigma_deriv_cl;

  results(check_y_not_nan, check_mu_finite, check_sigma_positive, lcdf_cl,
          y_deriv_cl, mu_deriv_cl, sigma_deriv_cl)
      = expressions(y_not_nan_expr, mu_finite_expr, sigma_positive_expr,
                    lcdf_expr, calc_if<!is_constant<T_y_cl>::value>(y_deriv),
                    calc_if<!is_constant<T_loc_cl>::value>(mu_deriv),
                    calc_if<!is_constant<T_scale_cl>::value>(sigma_deriv));

  double lcdf = sum(from_matrix_cl(lcdf_cl));

  auto ops_partials = make_partials_propagator(y_col, mu_col, sigma_col);

  if (!is_constant<T_y_cl>::value) {
    partials<0>(ops_partials) = std::move(y_deriv_cl);
  }
  if (!is_constant<T_loc_cl>::value) {
    partials<1>(ops_partials) = std::move(mu_deriv_cl);
  }
  if (!is_constant<T_scale_cl>::value) {
    partials<2>(ops_partials) = std::move(sigma_deriv_cl);
  }
  return ops_partials.build(lcdf);
}

}  // namespace math
}  // namespace stan
#endif
#endif
