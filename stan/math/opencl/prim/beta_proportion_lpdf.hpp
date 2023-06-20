#ifndef STAN_MATH_OPENCL_PRIM_BETA_PROPORTION_LPDF_HPP
#define STAN_MATH_OPENCL_PRIM_BETA_PROPORTION_LPDF_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1m.hpp>
#include <stan/math/prim/fun/log1p.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>

namespace stan {
namespace math {

/** \ingroup opencl
 * The log of the beta density for specified y, location, and
 * precision: beta_proportion_lpdf(y | mu, kappa) = beta_lpdf(y | mu *
 * kappa, (1 - mu) * kappa).  Any arguments other than scalars must be
 * containers of the same size.  With non-scalar arguments, the return
 * is the sum of the log pdfs with scalars broadcast as necessary.
 *
 * <p> The result log probability is defined to be the sum of
 * the log probabilities for each observation/mu/kappa triple.
 *
 * Prior location, mu, must be contained in (0, 1).  Prior precision
 * must be positive.
 *
 * @tparam T_y_cl type of outcome
 * @tparam T_loc_cl type of prior location
 * @tparam T_prec_cl type of prior precision
 *
 * @param y (Sequence of) dependant variable(s)
 * @param mu (Sequence of) location parameter(s)
 * @param kappa (Sequence of) precision parameter(s)
 * @return The log of the product of densities.
 */
template <bool propto, typename T_y_cl, typename T_loc_cl, typename T_prec_cl,
          require_all_prim_or_rev_kernel_expression_t<T_y_cl, T_loc_cl,
                                                      T_prec_cl>* = nullptr,
          require_any_not_stan_scalar_t<T_y_cl, T_loc_cl, T_prec_cl>* = nullptr>
return_type_t<T_y_cl, T_loc_cl, T_prec_cl> beta_proportion_lpdf(
    const T_y_cl& y, const T_loc_cl& mu, const T_prec_cl& kappa) {
  static const char* function = "beta_proportion_lpdf(OpenCL)";
  using T_partials_return = partials_return_t<T_y_cl, T_loc_cl, T_prec_cl>;
  using std::isfinite;

  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Precision parameter", kappa);
  const size_t N = max_size(y, mu, kappa);
  if (N == 0) {
    return 0.0;
  }
  if (!include_summand<propto, T_y_cl, T_loc_cl, T_prec_cl>::value) {
    return 0.0;
  }

  const auto& y_col = as_column_vector_or_scalar(y);
  const auto& mu_col = as_column_vector_or_scalar(mu);
  const auto& kappa_col = as_column_vector_or_scalar(kappa);

  const auto& y_val = value_of(y_col);
  const auto& mu_val = value_of(mu_col);
  const auto& kappa_val = value_of(kappa_col);

  auto check_y_bounded
      = check_cl(function, "Random variable", y_val, "in the interval [0, 1]");
  auto y_bounded_expr = 0 <= y_val && y_val <= 1;
  auto check_mu_bounded = check_cl(function, "Location parameter", mu_val,
                                   "in the interval (0, 1)");
  auto mu_bounded_expr = 0 < mu_val && mu_val < 1;
  auto check_kappa_positive_finite = check_cl(
      function, "Precision parameter", kappa_val, "in the interval [0, 1]");
  auto kappa_positive_finite = 0 < kappa_val && isfinite(kappa_val);

  auto log_y_expr = log(y_val);
  auto log1m_y_expr = log1p(-y_val);
  auto mukappa_expr = elt_multiply(mu_val, kappa_val);
  auto logp_expr = colwise_sum(
      elt_multiply(mukappa_expr - 1, log_y_expr)
      + elt_multiply(kappa_val - mukappa_expr - 1, log1m_y_expr)
      + static_select<include_summand<propto, T_prec_cl>::value>(
          lgamma(kappa_val), 0)
      - static_select<include_summand<propto, T_loc_cl, T_prec_cl>::value>(
          lgamma(mukappa_expr) + lgamma(kappa_val - mukappa_expr), 0));
  auto y_deriv_expr = elt_divide(mukappa_expr - 1, y_val)
                      + elt_divide(kappa_val - mukappa_expr - 1, y_val - 1);
  auto digamma_mukappa_expr = digamma(mukappa_expr);
  auto digamma_kappa_mukappa_expr = digamma(kappa_val - mukappa_expr);
  auto mu_deriv_expr = elt_multiply(kappa_val, digamma_kappa_mukappa_expr
                                                   - digamma_mukappa_expr
                                                   + log_y_expr - log1m_y_expr);
  auto kappa_deriv_expr
      = digamma(kappa_val)
        + elt_multiply(mu_val, log_y_expr - digamma_mukappa_expr)
        + elt_multiply(1 - mu_val, log1m_y_expr - digamma_kappa_mukappa_expr);

  matrix_cl<double> logp_cl;
  matrix_cl<double> y_deriv_cl;
  matrix_cl<double> mu_deriv_cl;
  matrix_cl<double> kappa_deriv_cl;

  results(check_y_bounded, check_mu_bounded, check_kappa_positive_finite,
          logp_cl, y_deriv_cl, mu_deriv_cl, kappa_deriv_cl)
      = expressions(y_bounded_expr, mu_bounded_expr, kappa_positive_finite,
                    logp_expr,
                    calc_if<!is_constant<T_y_cl>::value>(y_deriv_expr),
                    calc_if<!is_constant<T_loc_cl>::value>(mu_deriv_expr),
                    calc_if<!is_constant<T_prec_cl>::value>(kappa_deriv_expr));

  T_partials_return logp = sum(from_matrix_cl(logp_cl));

  auto ops_partials = make_partials_propagator(y_col, mu_col, kappa_col);
  if (!is_constant<T_y_cl>::value) {
    partials<0>(ops_partials) = std::move(y_deriv_cl);
  }
  if (!is_constant<T_loc_cl>::value) {
    partials<1>(ops_partials) = std::move(mu_deriv_cl);
  }
  if (!is_constant<T_prec_cl>::value) {
    partials<2>(ops_partials) = std::move(kappa_deriv_cl);
  }

  return ops_partials.build(logp);
}

}  // namespace math
}  // namespace stan
#endif
#endif
