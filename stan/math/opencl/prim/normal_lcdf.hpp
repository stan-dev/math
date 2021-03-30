#ifndef STAN_MATH_OPENCL_PRIM_NORMAL_LCDF_HPP
#define STAN_MATH_OPENCL_PRIM_NORMAL_LCDF_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/elt_divide.hpp>
#include <stan/math/prim/fun/elt_multiply.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>

namespace stan {
namespace math {

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
  static const char* function = "normal_lcdf(OpenCL)";
  using T_partials_return = partials_return_t<T_y_cl, T_loc_cl, T_scale_cl>;
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

  // for comments on numeric approximations used see CPU implementation
  auto scaled_diff = elt_divide(y_val - mu_val, sigma_val * SQRT_TWO);
  auto sigma_sqrt2 = sigma_val * SQRT_TWO;
  auto x2 = square(scaled_diff);

  auto x4 = pow(scaled_diff, 4);
  auto x6 = pow(scaled_diff, 6);
  auto x8 = pow(scaled_diff, 8);
  auto x10 = pow(scaled_diff, 10);
  auto temp_p = 0.000658749161529837803157
                + elt_divide(0.0160837851487422766278, x2)
                + elt_divide(0.125781726111229246204, x4)
                + elt_divide(0.360344899949804439429, x6)
                + elt_divide(0.305326634961232344035, x8)
                + elt_divide(0.0163153871373020978498, x10);
  auto temp_q = -0.00233520497626869185443
                - elt_divide(0.0605183413124413191178, x2)
                - elt_divide(0.527905102951428412248, x4)
                - elt_divide(1.87295284992346047209, x6)
                - elt_divide(2.56852019228982242072, x8) - elt_divide(1.0, x10);
  auto temp = LOG_HALF
              + log(INV_SQRT_PI + elt_divide(elt_divide(temp_p, temp_q), x2))
              - log(-scaled_diff) - x2;
  auto lcdf_n = select(
      scaled_diff > 0.0, log1p(-0.5 * erfc(scaled_diff)),
      select(scaled_diff > -20.0, log(erfc(-scaled_diff)) + LOG_HALF,
             select(10.0 * log(fabs(scaled_diff))
                        < log(std::numeric_limits<T_partials_return>::max()),
                    temp, NEGATIVE_INFTY)));
  auto lcdf_expr = colwise_sum(lcdf_n);

  auto c1 = scaled_diff > 2.9;
  auto c2 = scaled_diff > 2.5;
  auto c3 = scaled_diff > 2.1;
  auto c4 = scaled_diff > 1.5;
  auto c5 = scaled_diff > 0.8;
  auto c6 = scaled_diff > 0.1;
  auto c7 = 10.0 * log(fabs(scaled_diff))
            < log(std::numeric_limits<T_partials_return>::max());
  auto t = select(
      c1, elt_divide(1.0, 1.0 + 0.3275911 * scaled_diff),
      select(
          c2, scaled_diff - 2.7,
          select(c3, scaled_diff - 2.3,
                 select(c4, scaled_diff - 1.85,
                        select(c5, scaled_diff - 1.15,
                               select(c6, scaled_diff - 0.45,
                                      select(c7, elt_divide(
                                                 1.0,
                                                 1.0 - 0.3275911 * scaled_diff),
                                             0.0)))))));
  auto t2 = square(t);
  auto t4 = pow(t, 4);
  auto t6 = pow(t, 6);
  auto correction_term = select(
      scaled_diff < -29.0,
      0.0015065154280332 * x2 - 0.3993154819705530 * scaled_diff
          - 4.2919418242931700,
      select(scaled_diff < -17.0,
             0.0001263257217272 * elt_multiply(x2, scaled_diff)
                 + 0.0123586859488623 * x2 - 0.0860505264736028 * scaled_diff
                 - 1.252783383752970,
             select(scaled_diff < -7.0,
                    0.000471585349920831 * elt_multiply(x2, scaled_diff)
                        + 0.0296839305424034 * x2
                        + 0.207402143352332 * scaled_diff + 0.425316974683324,
                    select(scaled_diff < -3.9,
                           -0.0006972280656443 * elt_multiply(x2, scaled_diff)
                               + 0.0068218494628567 * x2
                               + 0.0585761964460277 * scaled_diff
                               + 0.1034397670201370, select(
                                   scaled_diff < -2.1,
                                   -0.0018742199480885
                                           * elt_multiply(x2, scaled_diff)
                                       - 0.0097119598291202 * x2
                                       - 0.0170137970924080 * scaled_diff
                                       - 0.0100428567412041,
                                   0.0)))));
  auto dncdf_log = select(
      c1,
      elt_divide(
          1.0, SQRT_PI
                   * (exp(x2) - 0.254829592 + 0.284496736 * t - 1.421413741 * t2
                      + 1.453152027 * elt_multiply(t2, t) - 1.061405429 * t4)),
      select(
          c2,
          0.0003849882382 - 0.002079084702 * t + 0.005229340880 * t2
              - 0.008029540137 * elt_multiply(t2, t) + 0.008232190507 * t4
              - 0.005692364250 * elt_multiply(t4, t) + 0.002399496363 * t6,
          select(
              c3,
              0.002846135439 - 0.01310032351 * t + 0.02732189391 * t2
                  - 0.03326906904 * elt_multiply(t2, t) + 0.02482478940 * t4
                  - 0.009883071924 * elt_multiply(t4, t) - 0.0002771362254 * t6,
              select(
                  c4,
                  0.01849212058 - 0.06876280470 * t + 0.1099906382 * t2
                      - 0.09274533184 * elt_multiply(t2, t) + 0.03543327418 * t4
                      + 0.005644855518 * elt_multiply(t4, t)
                      - 0.01111434424 * t6,
                  select(
                      c5,
                      0.1585747034 - 0.3898677543 * t + 0.3515963775 * t2
                          - 0.09748053605 * elt_multiply(t2, t)
                          - 0.04347986191 * t4
                          + 0.02182506378 * elt_multiply(t4, t)
                          + 0.01074751427 * t6,
                      select(
                          c6,
                          0.6245634904 - 0.9521866949 * t + 0.3986215682 * t2
                              + 0.04700850676 * elt_multiply(t2, t)
                              - 0.03478651979 * t4
                              - 0.01772675404 * elt_multiply(t4, t)
                              + 0.0006577254811 * t6,
                          select(
                              c7,
                              elt_divide(
                                  2.0,
                                  SQRT_PI
                                      * (0.254829592 * t - 0.284496736 * t2
                                         + 1.421413741 * elt_multiply(t2, t)
                                         - 1.453152027 * t4
                                         + 1.061405429 * elt_multiply(t4, t)))
                                  + correction_term,
                              INFTY)))))));
  auto y_deriv = elt_divide(dncdf_log, sigma_sqrt2);
  auto mu_deriv = -y_deriv;
  auto sigma_deriv
      = elt_divide(elt_multiply(dncdf_log, scaled_diff), sigma_val);

  matrix_cl<double> lcdf_cl;
  matrix_cl<double> y_deriv_cl;
  matrix_cl<double> mu_deriv_cl;
  matrix_cl<double> sigma_deriv_cl;

  results(check_y_not_nan, check_mu_finite, check_sigma_positive,
          lcdf_cl, y_deriv_cl, mu_deriv_cl, sigma_deriv_cl)
      = expressions(y_not_nan_expr, mu_finite_expr, sigma_positive_expr,
                    lcdf_expr,
                    calc_if<!is_constant<T_y_cl>::value>(y_deriv),
                    calc_if<!is_constant<T_loc_cl>::value>(mu_deriv),
                    calc_if<!is_constant<T_scale_cl>::value>(sigma_deriv));

  T_partials_return lcdf = sum(from_matrix_cl(lcdf_cl));

  operands_and_partials<decltype(y_col), decltype(mu_col), decltype(sigma_col)>
      ops_partials(y_col, mu_col, sigma_col);

  if (!is_constant<T_y_cl>::value) {
    ops_partials.edge1_.partials_ = std::move(y_deriv_cl);
  }
  if (!is_constant<T_loc_cl>::value) {
    ops_partials.edge2_.partials_ = std::move(mu_deriv_cl);
  }
  if (!is_constant<T_scale_cl>::value) {
    ops_partials.edge3_.partials_ = std::move(sigma_deriv_cl);
  }
  return ops_partials.build(lcdf);
}

}  // namespace math
}  // namespace stan
#endif
#endif
