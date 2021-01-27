#ifndef STAN_MATH_OPENCL_PRIM_BETA_BINOMIAL_LPMF_HPP
#define STAN_MATH_OPENCL_PRIM_BETA_BINOMIAL_LPMF_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/size.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>
#include <stan/math/prim/fun/elt_multiply.hpp>
#include <stan/math/prim/fun/elt_divide.hpp>

namespace stan {
namespace math {

/** \ingroup opencl
 * The log of the beta density for the specified scalar(s) given the specified
 * sample stan::math::size(s). y, alpha, or beta can each either be scalar or a
 * vector on OpenCL device. Any vector inputs must be the same length.
 *
 * <p> The result log probability is defined to be the sum of
 * the log probabilities for each observation/alpha/beta triple.
 *
 * Prior sample sizes, alpha and beta, must be greater than 0.
 *
 * @tparam T_n_cl type of scalar outcome
 * @tparam T_size1_cl type of prior scale for successes
 * @tparam T_size2_cl type of prior scale for failures
 * @param y (Sequence of) scalar(s).
 * @param alpha (Sequence of) prior sample stan::math::size(s).
 * @param beta (Sequence of) prior sample stan::math::size(s).
 * @return The log of the product of densities.
 */
template <
    bool propto, typename T_n_cl, typename T_N_cl, typename T_size1_cl,
    typename T_size2_cl,
    require_all_prim_or_rev_kernel_expression_t<T_n_cl, T_size1_cl,
                                                T_size2_cl>* = nullptr,
    require_any_not_stan_scalar_t<T_n_cl, T_size1_cl, T_size2_cl>* = nullptr>
return_type_t<T_n_cl, T_size1_cl, T_size2_cl> beta_binomial_lpmf(
    const T_n_cl& n, const T_N_cl N, const T_size1_cl& alpha,
    const T_size2_cl& beta) {
  using std::isfinite;
  static const char* function = "beta_binomial_lpmf(OpenCL)";

  check_consistent_sizes(function, "Successes variable", n,
                         "Population size parameter", N,
                         "First prior sample size parameter", alpha,
                         "Second prior sample size parameter", beta);
  const size_t N_size = max_size(n, N, alpha, beta);
  if (N_size == 0) {
    return 0.0;
  }
  if (!include_summand<propto, T_size1_cl, T_size2_cl>::value) {
    return 0.0;
  }

  const auto& alpha_val = value_of(alpha);
  const auto& beta_val = value_of(beta);

  auto check_N_nonnegative
      = check_cl(function, "Population size parameter", N, "nonnegative");
  auto N_nonnegative = N >= 0;
  auto check_alpha_pos_finite
      = check_cl(function, "First prior sample size parameter", alpha_val,
                 "positive finite");
  auto alpha_pos_finite = alpha_val > 0.0 && isfinite(alpha_val);
  auto check_beta_pos_finite
      = check_cl(function, "First prior sample size parameter", beta_val,
                 "positive finite");
  auto beta_pos_finite = beta_val > 0.0 && isfinite(beta_val);

  auto return_neg_inf = (n < 0 || n > N) + constant(0, N_size, 1);
  auto lbeta_diff
      = lbeta(n + alpha_val, N - n + beta_val) - lbeta(alpha_val, beta_val);
  auto digamma_diff
      = digamma(alpha_val + beta_val) - digamma(N + alpha_val + beta_val);
  auto logp_expr = colwise_sum(static_select<include_summand<propto>::value>(
      binomial_coefficient_log(N, n) + lbeta_diff, lbeta_diff));

  auto alpha_deriv = digamma(n + alpha_val) + digamma_diff - digamma(alpha_val);
  auto beta_deriv
      = digamma(N - n + beta_val) + digamma_diff - digamma(beta_val);

  matrix_cl<double> logp_cl;
  matrix_cl<int> return_neg_inf_cl;
  matrix_cl<double> alpha_deriv_cl;
  matrix_cl<double> beta_deriv_cl;

  results(check_N_nonnegative, check_alpha_pos_finite, check_beta_pos_finite,
          logp_cl, alpha_deriv_cl, beta_deriv_cl)
      = expressions(N_nonnegative, alpha_pos_finite, beta_pos_finite, logp_expr,
                    calc_if<!is_constant<T_size1_cl>::value>(alpha_deriv),
                    calc_if<!is_constant<T_size2_cl>::value>(beta_deriv));

  if (from_matrix_cl(return_neg_inf_cl).any()) {
    return LOG_ZERO;
  }

  double logp = sum(from_matrix_cl(logp_cl));

  operands_and_partials<T_size1_cl, T_size2_cl> ops_partials(alpha, beta);
  if (!is_constant<T_size1_cl>::value) {
    ops_partials.edge1_.partials_ = std::move(alpha_deriv_cl);
  }
  if (!is_constant<T_size2_cl>::value) {
    ops_partials.edge2_.partials_ = std::move(beta_deriv_cl);
  }

  return ops_partials.build(logp);
}

}  // namespace math
}  // namespace stan
#endif
#endif
