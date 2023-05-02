#ifndef STAN_MATH_OPENCL_PRIM_NEG_BINOMIAL_2_LPMF_HPP
#define STAN_MATH_OPENCL_PRIM_NEG_BINOMIAL_2_LPMF_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/elt_divide.hpp>
#include <stan/math/prim/fun/elt_multiply.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <stan/math/prim/fun/multiply_log.hpp>

namespace stan {
namespace math {

/** \ingroup opencl
 * The log of the negative binomial density for the specified scalars given
 * the specified mean(s) and deviation(s). n, mu, or phi can
 * each be either a scalar or a vector matrix_cl. Any vector inputs
 * must be the same length.
 *
 * <p>The result log probability is defined to be the sum of the
 * log probabilities for each observation/mean/deviation triple.
 *
 * @tparam T_n_cl type of scalar
 * @tparam T_location_cl type of location parameter
 * @tparam T_precision_cl type of precision parameter
 * @param n (Sequence of) scalar(s).
 * @param mu (Sequence of) location parameter(s)
 * @param phi (Sequence of) precision parameters
 * @return The log of the product of the densities.
 * @throw std::domain_error if the scale is not positive.
 */
template <bool propto, typename T_n_cl, typename T_location_cl,
          typename T_precision_cl,
          require_all_prim_or_rev_kernel_expression_t<
              T_n_cl, T_location_cl, T_precision_cl>* = nullptr,
          require_any_not_stan_scalar_t<T_n_cl, T_location_cl,
                                        T_precision_cl>* = nullptr>
inline return_type_t<T_n_cl, T_location_cl, T_precision_cl> neg_binomial_2_lpmf(
    const T_n_cl& n, const T_location_cl& mu, const T_precision_cl& phi) {
  static const char* function = "neg_binomial_2_lpmf(OpenCL)";
  using T_partials_return
      = partials_return_t<T_n_cl, T_location_cl, T_precision_cl>;
  using std::isfinite;
  using std::isnan;

  check_consistent_sizes(function, "Failures variable", n, "Location parameter",
                         mu, "Precision parameter", phi);
  const size_t N = max_size(n, mu, phi);
  if (N == 0) {
    return 0.0;
  }
  if (!include_summand<propto, T_n_cl, T_location_cl, T_precision_cl>::value) {
    return 0.0;
  }

  const auto& mu_col = as_column_vector_or_scalar(mu);
  const auto& phi_col = as_column_vector_or_scalar(phi);

  const auto& mu_val = value_of(mu_col);
  const auto& phi_val = value_of(phi_col);

  auto check_n_nonnegative
      = check_cl(function, "Failures variable", n, "nonnegative");
  auto n_nonnegative = n >= 0;
  auto check_mu_positive_finite
      = check_cl(function, "Log location parameter", mu_val, "positive finite");
  auto mu_positive_finite = 0 < mu_val && isfinite(mu_val);
  auto check_phi_positive_finite
      = check_cl(function, "Precision parameter", phi_val, "positive finite");
  auto phi_positive_finite = 0 < phi_val && isfinite(phi_val);

  auto log_phi = log(phi_val);
  auto mu_plus_phi = mu_val + phi_val;
  auto log_mu_plus_phi = log(mu_plus_phi);
  auto n_plus_phi = n + phi_val;

  auto logp1 = -elt_multiply(phi_val, log1p(elt_divide(mu_val, phi_val)))
               - elt_multiply(n, log_mu_plus_phi);
  auto logp2 = static_select<include_summand<propto, T_precision_cl>::value>(
      logp1 + binomial_coefficient_log(n_plus_phi - 1, n), logp1);
  auto logp_expr = colwise_sum(
      static_select<include_summand<propto, T_location_cl>::value>(
          logp2 + multiply_log(n, mu_val), logp2));

  auto mu_deriv = elt_divide(n, mu_val) - elt_divide(n + phi_val, mu_plus_phi);
  auto log_term
      = select(mu_val < phi_val, log1p(-elt_divide(mu_val, mu_plus_phi)),
               log_phi - log_mu_plus_phi);
  auto phi_deriv = elt_divide(mu_val - n, mu_plus_phi) + log_term
                   - digamma(phi_val) + digamma(n_plus_phi);

  matrix_cl<double> logp_cl;
  matrix_cl<double> mu_deriv_cl;
  matrix_cl<double> phi_deriv_cl;

  results(check_n_nonnegative, check_mu_positive_finite,
          check_phi_positive_finite, logp_cl, mu_deriv_cl, phi_deriv_cl)
      = expressions(n_nonnegative, mu_positive_finite, phi_positive_finite,
                    logp_expr,
                    calc_if<!is_constant<T_location_cl>::value>(mu_deriv),
                    calc_if<!is_constant<T_precision_cl>::value>(phi_deriv));

  T_partials_return logp = sum(from_matrix_cl(logp_cl));

  auto ops_partials = make_partials_propagator(mu_col, phi_col);

  if (!is_constant<T_location_cl>::value) {
    partials<0>(ops_partials) = std::move(mu_deriv_cl);
  }
  if (!is_constant<T_precision_cl>::value) {
    partials<1>(ops_partials) = std::move(phi_deriv_cl);
  }
  return ops_partials.build(logp);
}

}  // namespace math
}  // namespace stan
#endif
#endif
