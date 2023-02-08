#ifndef STAN_MATH_OPENCL_PRIM_POISSON_LPMF_HPP
#define STAN_MATH_OPENCL_PRIM_POISSON_LPMF_HPP
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
 * Returns the log PMF of the Poisson distribution. If containers are
 * supplied, returns the log sum of the probabilities.
 *
 * @tparam T_n_cl type of integer parameters
 * @tparam T_rate_cl type of chance of success parameters
 * @param n integer parameter
 * @param lambda rate parameter
 * @return log probability or log sum of probabilities
 * @throw std::domain_error if lambda is not a valid probability
 * @throw std::invalid_argument if container sizes mismatch.
 */
template <
    bool propto, typename T_n_cl, typename T_rate_cl,
    require_all_prim_or_rev_kernel_expression_t<T_n_cl, T_rate_cl>* = nullptr,
    require_any_not_stan_scalar_t<T_n_cl, T_rate_cl>* = nullptr>
return_type_t<T_rate_cl> poisson_lpmf(const T_n_cl& n,
                                      const T_rate_cl& lambda) {
  static const char* function = "poisson_lpmf(OpenCL)";
  using T_partials_return = partials_return_t<T_rate_cl>;
  using std::isinf;
  constexpr bool is_n_vector = !is_stan_scalar<T_n_cl>::value;
  constexpr bool is_lambda_vector = !is_stan_scalar<T_rate_cl>::value;

  check_consistent_sizes(function, "Random variable", n, "Rate parameter",
                         lambda);
  const size_t N = is_n_vector ? math::size(n) : math::size(lambda);
  if (N == 0) {
    return 0.0;
  }
  if (!include_summand<propto, T_rate_cl>::value) {
    return 0.0;
  }

  const auto& lambda_col = as_column_vector_or_scalar(lambda);
  const auto& lambda_val = value_of(lambda_col);

  T_partials_return logp(0.0);
  operands_and_partials<decltype(lambda_col)> ops_partials(lambda_col);

  auto check_n_nonnegative
      = check_cl(function, "Random variable", n, "nonnegative");
  auto n_nonnegative = 0 <= n;
  auto check_lambda_nonnegative
      = check_cl(function, "Log rate parameter", lambda_val, "nonnegative");
  auto lambda_nonnegative = 0.0 <= lambda_val;

  auto return_log_zero = colwise_max(
      cast<char>(isinf(lambda_val) || (lambda_val == 0 && n != 0)));

  auto logp1 = multiply_log(n, lambda_val);
  auto logp2 = static_select<include_summand<propto, T_rate_cl>::value>(
      logp1 - lambda_val, logp1);
  auto logp_expr = colwise_sum(static_select<include_summand<propto>::value>(
      logp2 - lgamma(n + 1.0), logp2));

  auto deriv = elt_divide(n, lambda_val) - 1.0;

  matrix_cl<char> return_log_zero_cl;
  matrix_cl<double> logp_cl;
  matrix_cl<double> deriv_cl;

  results(check_n_nonnegative, check_lambda_nonnegative, return_log_zero_cl,
          logp_cl, deriv_cl)
      = expressions(n_nonnegative, lambda_nonnegative, return_log_zero,
                    logp_expr,
                    calc_if<!is_constant_all<T_rate_cl>::value>(deriv));

  if (from_matrix_cl(return_log_zero_cl).any()) {
    return LOG_ZERO;
  }

  logp = sum(from_matrix_cl(logp_cl));

  if (!is_constant_all<T_rate_cl>::value) {
    ops_partials.edge1_.partials_ = deriv_cl;
  }

  return ops_partials.build(logp);
}

}  // namespace math
}  // namespace stan
#endif
#endif
