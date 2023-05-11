#ifndef STAN_MATH_OPENCL_PRIM_STD_NORMAL_LCDF_HPP
#define STAN_MATH_OPENCL_PRIM_STD_NORMAL_LCDF_HPP
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
const char opencl_std_normal_lcdf_impl[] = STRINGIFY(
    double std_normal_lcdf_lcdf_n; if (std_normal_lcdf_scaled_y > 0.0) {
      // CDF(x) = 1/2 + 1/2erf(x) = 1 - 1/2erfc(x)
      std_normal_lcdf_lcdf_n = log1p(-0.5 * erfc(std_normal_lcdf_scaled_y));
      if (isnan(std_normal_lcdf_lcdf_n)) {
        std_normal_lcdf_lcdf_n = 0;
      }
    } else if (std_normal_lcdf_scaled_y > -20.0) {
      // CDF(x) = 1/2 - 1/2erf(-x) = 1/2erfc(-x)
      std_normal_lcdf_lcdf_n = log(erfc(-std_normal_lcdf_scaled_y)) - M_LN2;
    } else if (10.0 * log(fabs(std_normal_lcdf_scaled_y)) < log(DBL_MAX)) {
      // entering territory where erfc(-x)~0
      // need to use direct numerical approximation of lcdf instead
      // the following based on W. J. Cody, Math. Comp. 23(107):631-638 (1969)
      // CDF(x) = 1/2erfc(-x)
      const double x4 = pow(std_normal_lcdf_scaled_y, 4);
      const double x6 = pow(std_normal_lcdf_scaled_y, 6);
      const double x8 = pow(std_normal_lcdf_scaled_y, 8);
      const double x10 = pow(std_normal_lcdf_scaled_y, 10);
      const double temp_p
          = 0.000658749161529837803157
            + 0.0160837851487422766278 / std_normal_lcdf_x2
            + 0.125781726111229246204 / x4 + 0.360344899949804439429 / x6
            + 0.305326634961232344035 / x8 + 0.0163153871373020978498 / x10;
      const double temp_q = -0.00233520497626869185443
                            - 0.0605183413124413191178 / std_normal_lcdf_x2
                            - 0.527905102951428412248 / x4
                            - 1.87295284992346047209 / x6
                            - 2.56852019228982242072 / x8 - 1.0 / x10;
      std_normal_lcdf_lcdf_n
          += log(0.5 * M_2_SQRTPI + (temp_p / temp_q) / std_normal_lcdf_x2)
             - M_LN2 - log(-std_normal_lcdf_scaled_y) - std_normal_lcdf_x2;
    } else {
      // std_normal_lcdf_scaled_y^10 term will overflow
      std_normal_lcdf_lcdf_n = -INFINITY;
    });
const char opencl_std_normal_lcdf_dnlcdf[] = STRINGIFY(
    // compute partial derivatives
    // based on analytic form given by:
    // dln(CDF)/dx = exp(-x^2)/(sqrt(pi)*(1/2+erf(x)/2)
    double std_normal_lcdf_dnlcdf = 0.0; double t = 0.0; double t2 = 0.0;
    double t4 = 0.0;

    // calculate using piecewise function
    // (due to instability / inaccuracy in the various approximations)
    if (std_normal_lcdf_deriv_scaled_y > 2.9) {
      // approximation derived from Abramowitz and Stegun (1964) 7.1.26
      t = 1.0 / (1.0 + 0.3275911 * std_normal_lcdf_deriv_scaled_y);
      t2 = t * t;
      t4 = pow(t, 4);
      std_normal_lcdf_dnlcdf
          = 0.5 * M_2_SQRTPI
            / (exp(std_normal_lcdf_deriv_x2) - 0.254829592 + 0.284496736 * t
               - 1.421413741 * t2 + 1.453152027 * t2 * t - 1.061405429 * t4);
    } else if (std_normal_lcdf_deriv_scaled_y > 2.5) {
      // in the trouble area where all of the standard numerical
      // approximations are unstable - bridge the gap using Taylor
      // expansions of the analytic function
      // use Taylor expansion centred around x=2.7
      t = std_normal_lcdf_deriv_scaled_y - 2.7;
      t2 = t * t;
      t4 = pow(t, 4);
      std_normal_lcdf_dnlcdf = 0.0003849882382 - 0.002079084702 * t
                               + 0.005229340880 * t2 - 0.008029540137 * t2 * t
                               + 0.008232190507 * t4 - 0.005692364250 * t4 * t
                               + 0.002399496363 * pow(t, 6);
    } else if (std_normal_lcdf_deriv_scaled_y > 2.1) {
      // use Taylor expansion centred around x=2.3
      t = std_normal_lcdf_deriv_scaled_y - 2.3;
      t2 = t * t;
      t4 = pow(t, 4);
      std_normal_lcdf_dnlcdf = 0.002846135439 - 0.01310032351 * t
                               + 0.02732189391 * t2 - 0.03326906904 * t2 * t
                               + 0.02482478940 * t4 - 0.009883071924 * t4 * t
                               - 0.0002771362254 * pow(t, 6);
    } else if (std_normal_lcdf_deriv_scaled_y > 1.5) {
      // use Taylor expansion centred around x=1.85
      t = std_normal_lcdf_deriv_scaled_y - 1.85;
      t2 = t * t;
      t4 = pow(t, 4);
      std_normal_lcdf_dnlcdf = 0.01849212058 - 0.06876280470 * t
                               + 0.1099906382 * t2 - 0.09274533184 * t2 * t
                               + 0.03543327418 * t4 + 0.005644855518 * t4 * t
                               - 0.01111434424 * pow(t, 6);
    } else if (std_normal_lcdf_deriv_scaled_y > 0.8) {
      // use Taylor expansion centred around x=1.15
      t = std_normal_lcdf_deriv_scaled_y - 1.15;
      t2 = t * t;
      t4 = pow(t, 4);
      std_normal_lcdf_dnlcdf = 0.1585747034 - 0.3898677543 * t
                               + 0.3515963775 * t2 - 0.09748053605 * t2 * t
                               - 0.04347986191 * t4 + 0.02182506378 * t4 * t
                               + 0.01074751427 * pow(t, 6);
    } else if (std_normal_lcdf_deriv_scaled_y > 0.1) {
      // use Taylor expansion centred around x=0.45
      t = std_normal_lcdf_deriv_scaled_y - 0.45;
      t2 = t * t;
      t4 = pow(t, 4);
      std_normal_lcdf_dnlcdf = 0.6245634904 - 0.9521866949 * t
                               + 0.3986215682 * t2 + 0.04700850676 * t2 * t
                               - 0.03478651979 * t4 - 0.01772675404 * t4 * t
                               + 0.0006577254811 * pow(t, 6);
    } else if (10.0 * log(fabs(std_normal_lcdf_deriv_scaled_y))
               < log(DBL_MAX)) {
      // approximation derived from Abramowitz and Stegun (1964) 7.1.26
      // use fact that erf(x)=-erf(-x)
      // Abramowitz and Stegun define this for -inf<x<0 but seems to be
      // accurate for -inf<x<0.1
      t = 1.0 / (1.0 - 0.3275911 * std_normal_lcdf_deriv_scaled_y);
      t2 = t * t;
      t4 = pow(t, 4);
      std_normal_lcdf_dnlcdf
          = M_2_SQRTPI
            / (0.254829592 * t - 0.284496736 * t2 + 1.421413741 * t2 * t
               - 1.453152027 * t4 + 1.061405429 * t4 * t);
      // check if we need to add a correction term
      // (from cubic fit of residuals)
      if (std_normal_lcdf_deriv_scaled_y < -29.0) {
        std_normal_lcdf_dnlcdf
            += 0.0015065154280332 * std_normal_lcdf_deriv_x2
               - 0.3993154819705530 * std_normal_lcdf_deriv_scaled_y
               - 4.2919418242931700;
      } else if (std_normal_lcdf_deriv_scaled_y < -17.0) {
        std_normal_lcdf_dnlcdf
            += 0.0001263257217272 * std_normal_lcdf_deriv_x2
                   * std_normal_lcdf_deriv_scaled_y
               + 0.0123586859488623 * std_normal_lcdf_deriv_x2
               - 0.0860505264736028 * std_normal_lcdf_deriv_scaled_y
               - 1.252783383752970;
      } else if (std_normal_lcdf_deriv_scaled_y < -7.0) {
        std_normal_lcdf_dnlcdf
            += 0.000471585349920831 * std_normal_lcdf_deriv_x2
                   * std_normal_lcdf_deriv_scaled_y
               + 0.0296839305424034 * std_normal_lcdf_deriv_x2
               + 0.207402143352332 * std_normal_lcdf_deriv_scaled_y
               + 0.425316974683324;
      } else if (std_normal_lcdf_deriv_scaled_y < -3.9) {
        std_normal_lcdf_dnlcdf
            += -0.0006972280656443 * std_normal_lcdf_deriv_x2
                   * std_normal_lcdf_deriv_scaled_y
               + 0.0068218494628567 * std_normal_lcdf_deriv_x2
               + 0.0585761964460277 * std_normal_lcdf_deriv_scaled_y
               + 0.1034397670201370;
      } else if (std_normal_lcdf_deriv_scaled_y < -2.1) {
        std_normal_lcdf_dnlcdf
            += -0.0018742199480885 * std_normal_lcdf_deriv_x2
                   * std_normal_lcdf_deriv_scaled_y
               - 0.0097119598291202 * std_normal_lcdf_deriv_x2
               - 0.0170137970924080 * std_normal_lcdf_deriv_scaled_y
               - 0.0100428567412041;
      }
    } else { std_normal_lcdf_dnlcdf = INFINITY; });
}  // namespace internal

/** \ingroup opencl
 * Returns the log standard normal complementary cumulative distribution
 * function.
 *
 * @tparam T_y_cl type of scalar outcome
 * @param y (Sequence of) scalar(s).
 * @return The log of the product of densities.
 */
template <typename T_y_cl,
          require_all_prim_or_rev_kernel_expression_t<T_y_cl>* = nullptr,
          require_any_not_stan_scalar_t<T_y_cl>* = nullptr>
return_type_t<T_y_cl> std_normal_lcdf(const T_y_cl& y) {
  static const char* function = "std_normal_lcdf(OpenCL)";
  using std::isfinite;
  using std::isnan;

  const size_t N = math::size(y);
  if (N == 0) {
    return 1.0;
  }

  const auto& y_col = as_column_vector_or_scalar(y);
  const auto& y_val = value_of(y_col);

  auto check_y_not_nan
      = check_cl(function, "Random variable", y_val, "not NaN");
  auto y_not_nan_expr = !isnan(y_val);

  auto scaled_y = y_val * INV_SQRT_TWO;
  auto x2 = square(scaled_y);
  auto lcdf_expr = colwise_sum(
      opencl_code<internal::opencl_std_normal_lcdf_impl>(
          std::make_tuple("std_normal_lcdf_scaled_y", "std_normal_lcdf_x2"),
          scaled_y, x2)
          .template output<double>("std_normal_lcdf_lcdf_n"));
  auto dnlcdf = opencl_code<internal::opencl_std_normal_lcdf_dnlcdf>(
                    std::make_tuple("std_normal_lcdf_deriv_scaled_y",
                                    "std_normal_lcdf_deriv_x2"),
                    scaled_y, x2)
                    .template output<double>("std_normal_lcdf_dnlcdf");
  auto y_deriv = dnlcdf * INV_SQRT_TWO;

  matrix_cl<double> lcdf_cl;
  matrix_cl<double> y_deriv_cl;

  results(check_y_not_nan, lcdf_cl, y_deriv_cl) = expressions(
      y_not_nan_expr, lcdf_expr, calc_if<!is_constant<T_y_cl>::value>(y_deriv));

  double lcdf = from_matrix_cl(lcdf_cl).sum();

  auto ops_partials = make_partials_propagator(y_col);

  if (!is_constant<T_y_cl>::value) {
    partials<0>(ops_partials) = std::move(y_deriv_cl);
  }
  return ops_partials.build(lcdf);
}

}  // namespace math
}  // namespace stan
#endif
#endif
