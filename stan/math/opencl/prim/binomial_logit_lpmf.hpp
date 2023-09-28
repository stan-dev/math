#ifndef STAN_MATH_OPENCL_PRIM_BINOMIAL_LOGIT_LPMF_HPP
#define STAN_MATH_OPENCL_PRIM_BINOMIAL_LOGIT_LPMF_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <stan/math/prim/fun/binomial_coefficient_log.hpp>
#include <stan/math/prim/fun/log_inv_logit.hpp>
#include <stan/math/prim/fun/log1m_inv_logit.hpp>

namespace stan {
namespace math {

/** \ingroup opencl
 * Binomial log PMF in logit parametrization. Binomial(n|n, inv_logit(alpha))
 *
 * If given vectors of matching lengths, returns
 * the log sum of probabilities.
 *
 * @param n successes variable
 * @param N population size parameter
 * @param alpha logit transformed probability parameter
 * @return log probability or log sum of probabilities
 * @throw std::domain_error if N is negative or probability parameter is invalid
 * @throw std::invalid_argument if vector sizes do not match
 */
template <bool propto, typename T_n_cl, typename T_N_cl, typename T_prob_cl,
          require_all_prim_or_rev_kernel_expression_t<T_n_cl, T_N_cl,
                                                      T_prob_cl>* = nullptr,
          require_any_nonscalar_prim_or_rev_kernel_expression_t<
              T_n_cl, T_N_cl, T_prob_cl>* = nullptr>
return_type_t<T_prob_cl> binomial_logit_lpmf(const T_n_cl& n, const T_N_cl N,
                                             const T_prob_cl& alpha) {
  static const char* function = "binomial_logit_lpmf(OpenCL)";
  using T_partials_return = partials_return_t<T_prob_cl>;
  using std::isfinite;

  check_consistent_sizes(function, "Successes variable", n,
                         "Population size parameter", N,
                         "Probability parameter", alpha);
  const size_t siz = max_size(n, N, alpha);
  if (siz == 0) {
    return 0.0;
  }
  if (!include_summand<propto, T_prob_cl>::value) {
    return 0.0;
  }

  const auto& alpha_col = as_column_vector_or_scalar(alpha);
  const auto& alpha_val = value_of(alpha_col);

  auto check_n_bounded
      = check_cl(function, "Successes variable", n, "in the interval [0, N]");
  auto n_bounded = 0 <= n && n <= N;
  auto check_N_nonnegative
      = check_cl(function, "Population size variable", n, "nonnegative");
  auto N_nonnegative = N >= 0;
  auto check_alpha_finite
      = check_cl(function, "Probability parameter", alpha_val, "finite");
  auto alpha_finite = isfinite(alpha_val);

  auto log_inv_logit_alpha = log_inv_logit(alpha_val);
  auto log1m_inv_logit_alpha = log1m_inv_logit(alpha_val);
  auto n_diff = N - n;
  auto logp_expr1 = elt_multiply(n, log_inv_logit_alpha)
                    + elt_multiply(n_diff, log1m_inv_logit_alpha);
  auto logp_expr
      = static_select<include_summand<propto, T_n_cl, T_N_cl>::value>(
          logp_expr1 + binomial_coefficient_log(N, n), logp_expr1);
  auto alpha_deriv = n - elt_multiply(N, exp(log_inv_logit_alpha));

  matrix_cl<double> logp_cl;
  matrix_cl<double> alpha_deriv_cl;

  results(check_n_bounded, check_N_nonnegative, check_alpha_finite, logp_cl,
          alpha_deriv_cl)
      = expressions(n_bounded, N_nonnegative, alpha_finite, logp_expr,
                    calc_if<!is_constant<T_prob_cl>::value>(alpha_deriv));

  T_partials_return logp = sum(from_matrix_cl(logp_cl));
  auto ops_partials = make_partials_propagator(alpha_col);
  if (!is_constant<T_prob_cl>::value) {
    partials<0>(ops_partials) = std::move(alpha_deriv_cl);
  }

  return ops_partials.build(logp);
}

}  // namespace math
}  // namespace stan
#endif
#endif
