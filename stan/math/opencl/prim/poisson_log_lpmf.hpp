#ifndef STAN_MATH_OPENCL_PRIM_POISSON_LOG_LPMF_HPP
#define STAN_MATH_OPENCL_PRIM_POISSON_LOG_LPMF_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>

namespace stan {
namespace math {

/** \ingroup opencl
 * Returns the log PMF of the Poisson log distribution. If containers are
 * supplied, returns the log sum of the probabilities.
 *
 * @tparam T_n_cl type of integer parameters
 * @tparam T_log_rate_cl type of chance of success parameters
 * @param n integer parameter
 * @param alpha log rate parameter
 * @return log probability or log sum of probabilities
 * @throw std::domain_error if alpha is not a valid probability
 * @throw std::invalid_argument if container sizes mismatch.
 */
template <bool propto, typename T_n_cl, typename T_log_rate_cl,
          require_all_prim_or_rev_kernel_expression_t<T_n_cl,
                                                      T_log_rate_cl>* = nullptr,
          require_any_not_stan_scalar_t<T_n_cl, T_log_rate_cl>* = nullptr>
return_type_t<T_log_rate_cl> poisson_log_lpmf(const T_n_cl& n,
                                              const T_log_rate_cl& alpha) {
  static const char* function = "poisson_log_lpmf(OpenCL)";
  using T_partials_return = partials_return_t<T_log_rate_cl>;
  using std::isinf;
  using std::isnan;
  constexpr bool is_n_vector = !is_stan_scalar<T_n_cl>::value;
  constexpr bool is_alpha_vector = !is_stan_scalar<T_log_rate_cl>::value;

  check_consistent_sizes(function, "Random variable", n, "Log rate parameter",
                         alpha);
  const size_t N = is_n_vector ? math::size(n) : math::size(alpha);
  if (N == 0) {
    return 0.0;
  }
  if (!include_summand<propto, T_log_rate_cl>::value) {
    return 0.0;
  }

  const auto& alpha_col = as_column_vector_or_scalar(alpha);
  const auto& alpha_val = value_of(alpha_col);

  T_partials_return logp(0.0);
  operands_and_partials<decltype(alpha_col)> ops_partials(alpha_col);

  auto check_n_nonnegative
      = check_cl(function, "Random variable", n, "nonnegative");
  auto n_nonnegativer = 0 <= n;
  auto check_alpha_not_nan
      = check_cl(function, "Log rate parameter", alpha_val, "not nan");
  auto alpha_not_nan = !isnan(alpha_val);

  auto return_log_zero
      = colwise_max(cast<char>(isinf(alpha_val) && (alpha_val > 0 || n != 0)));
  auto exp_alpha = exp(alpha_val);

  auto logp1 = elt_multiply(n, alpha_val);
  auto logp2 = static_select<include_summand<propto, T_log_rate_cl>::value>(
      logp1 - exp_alpha, logp1);
  auto logp_expr = colwise_sum(static_select<include_summand<propto>::value>(
      logp2 - lgamma(n + 1.0), logp2));

  auto deriv = n - exp_alpha;

  matrix_cl<char> return_log_zero_cl;
  matrix_cl<double> logp_cl;
  matrix_cl<double> deriv_cl;

  results(check_n_nonnegative, check_alpha_not_nan, return_log_zero_cl, logp_cl,
          deriv_cl)
      = expressions(n_nonnegativer, alpha_not_nan, return_log_zero, logp_expr,
                    calc_if<!is_constant_all<T_log_rate_cl>::value>(deriv));

  if (from_matrix_cl(return_log_zero_cl).any()) {
    return LOG_ZERO;
  }

  logp = sum(from_matrix_cl(logp_cl));

  if (!is_constant_all<T_log_rate_cl>::value) {
    ops_partials.edge1_.partials_ = deriv_cl;
  }

  return ops_partials.build(logp);
}

}  // namespace math
}  // namespace stan
#endif
#endif
