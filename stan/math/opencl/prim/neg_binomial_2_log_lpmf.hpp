#ifndef STAN_MATH_OPENCL_PRIM_NEG_BINOMIAL_2_LOG_LPMF_HPP
#define STAN_MATH_OPENCL_PRIM_NEG_BINOMIAL_2_LOG_LPMF_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/elt_divide.hpp>
#include <stan/math/prim/fun/elt_multiply.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <stan/math/prim/fun/log1p_exp.hpp>

namespace stan {
namespace math {

/** \ingroup opencl
 * The log of the log transformed negative binomial density for the specified
 * scalars given the specified mean(s) and deviation(s). n, eta, or phi can each
 * be either a scalar or a vector matrix_cl. Any vector inputs must be the same
 * length.
 *
 * <p>The result log probability is defined to be the sum of the
 * log probabilities for each observation/mean/deviation triple.
 *
 * @tparam T_n_cl type of scalar
 * @tparam T_log_location_cl type of location parameter
 * @tparam T_precision_cl type of precision parameter
 * @param n (Sequence of) scalar(s).
 * @param eta (Sequence of) location parameter(s)
 * @param phi (Sequence of) precision parameters
 * @return The log of the product of the densities.
 * @throw std::domain_error if the scale is not positive.
 */
template <bool propto, typename T_n_cl, typename T_log_location_cl,
          typename T_precision_cl,
          require_all_prim_or_rev_kernel_expression_t<
              T_n_cl, T_log_location_cl, T_precision_cl>* = nullptr,
          require_any_not_stan_scalar_t<T_n_cl, T_log_location_cl,
                                        T_precision_cl>* = nullptr>
inline return_type_t<T_n_cl, T_log_location_cl, T_precision_cl>
neg_binomial_2_log_lpmf(const T_n_cl& n, const T_log_location_cl& eta,
                        const T_precision_cl& phi) {
  static const char* function = "neg_binomial_2_log_lpmf(OpenCL)";
  using T_partials_return
      = partials_return_t<T_n_cl, T_log_location_cl, T_precision_cl>;
  using std::isfinite;
  using std::isnan;

  check_consistent_sizes(function, "Failures variable", n,
                         "Log location parameter", eta, "Precision parameter",
                         phi);
  const size_t N = max_size(n, eta, phi);
  if (N == 0) {
    return 0.0;
  }
  if (!include_summand<propto, T_n_cl, T_log_location_cl,
                       T_precision_cl>::value) {
    return 0.0;
  }

  const auto& eta_col = as_column_vector_or_scalar(eta);
  const auto& phi_col = as_column_vector_or_scalar(phi);

  const auto& eta_val = value_of(eta_col);
  const auto& phi_val = value_of(phi_col);

  auto check_n_nonnegative
      = check_cl(function, "Failures variable", n, "nonnegative");
  auto n_nonnegative = n >= 0;
  auto check_eta_finite
      = check_cl(function, "Log location parameter", eta_val, "finite");
  auto eta_finite = isfinite(eta_val);
  auto check_phi_positive_finite
      = check_cl(function, "Precision parameter", phi_val, "positive finite");
  auto phi_positive_finite = 0 < phi_val && isfinite(phi_val);

  auto log_phi = log(phi_val);
  auto exp_eta = exp(eta_val);
  auto exp_eta_over_exp_eta_phi
      = elt_divide(1.0, elt_divide(phi_val, exp_eta) + 1.0);
  auto log1p_exp_eta_m_logphi = log1p_exp(eta_val - log_phi);
  auto n_plus_phi = n + phi_val;

  auto logp1 = -elt_multiply(phi_val, log1p_exp_eta_m_logphi)
               - elt_multiply(n, log_phi + log1p_exp_eta_m_logphi);
  auto logp2 = static_select<include_summand<propto, T_precision_cl>::value>(
      logp1 + binomial_coefficient_log(n_plus_phi - 1, n), logp1);
  auto logp_expr = colwise_sum(
      static_select<include_summand<propto, T_log_location_cl>::value>(
          logp2 + elt_multiply(n, eta_val), logp2));

  auto eta_deriv = n - elt_multiply(n_plus_phi, exp_eta_over_exp_eta_phi);
  auto phi_deriv = exp_eta_over_exp_eta_phi - elt_divide(n, exp_eta + phi_val)
                   - log1p_exp_eta_m_logphi - digamma(phi_val)
                   + digamma(n_plus_phi);

  matrix_cl<double> logp_cl;
  matrix_cl<double> eta_deriv_cl;
  matrix_cl<double> phi_deriv_cl;

  results(check_n_nonnegative, check_eta_finite, check_phi_positive_finite,
          logp_cl, eta_deriv_cl, phi_deriv_cl)
      = expressions(n_nonnegative, eta_finite, phi_positive_finite, logp_expr,
                    calc_if<!is_constant<T_log_location_cl>::value>(eta_deriv),
                    calc_if<!is_constant<T_precision_cl>::value>(phi_deriv));

  T_partials_return logp = sum(from_matrix_cl(logp_cl));

  auto ops_partials = make_partials_propagator(eta_col, phi_col);

  if (!is_constant<T_log_location_cl>::value) {
    partials<0>(ops_partials) = std::move(eta_deriv_cl);
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
