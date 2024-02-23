#ifndef STAN_MATH_OPENCL_PRIM_BINOMIAL_LPMF_HPP
#define STAN_MATH_OPENCL_PRIM_BINOMIAL_LPMF_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <stan/math/prim/fun/binomial_coefficient_log.hpp>

namespace stan {
namespace math {

/** \ingroup opencl
 * Returns the log PMF for the binomial distribution evaluated at the
 * specified success, population size, and chance of success. If given
 * containers of matching lengths, returns the log sum of probabilities.
 *
 * @tparam T_n_cl type of successes parameter
 * @tparam T_N_cl type of population size parameter
 * @tparam T_prob_cl type of chance of success parameter
 * @param n successes parameter
 * @param N population size parameter
 * @param theta chance of success parameter
 * @return log sum of probabilities
 * @throw std::domain_error if n is negative or greater than N
 * @throw std::domain_error if N is negative
 * @throw std::domain_error if theta is not a valid probability
 * @throw std::invalid_argument if container sizes mismatch
 */
template <bool propto, typename T_n_cl, typename T_N_cl, typename T_prob_cl,
          require_all_prim_or_rev_kernel_expression_t<T_n_cl, T_N_cl,
                                                      T_prob_cl>* = nullptr,
          require_any_nonscalar_prim_or_rev_kernel_expression_t<
              T_n_cl, T_N_cl, T_prob_cl>* = nullptr>
return_type_t<T_prob_cl> binomial_lpmf(const T_n_cl& n, const T_N_cl N,
                                       const T_prob_cl& theta) {
  static constexpr const char* function = "binomial_lpmf(OpenCL)";
  using T_partials_return = partials_return_t<T_prob_cl>;
  using std::isnan;

  check_consistent_sizes(function, "Successes variable", n,
                         "Population size parameter", N,
                         "Probability parameter", theta);
  const size_t siz = max_size(n, N, theta);
  if (siz == 0) {
    return 0.0;
  }
  if (!include_summand<propto, T_prob_cl>::value) {
    return 0.0;
  }

  const auto& theta_col = as_column_vector_or_scalar(theta);
  const auto& theta_val = value_of(theta_col);

  auto check_n_bounded
      = check_cl(function, "Successes variable", n, "in the interval [0, N]");
  auto n_bounded = 0 <= n && n <= N;
  auto check_N_nonnegative
      = check_cl(function, "Population size variable", n, "nonnegative");
  auto N_nonnegative = N >= 0;
  auto check_theta_bounded = check_cl(function, "Probability parameter",
                                      theta_val, "in the interval [0, 1]");
  auto theta_bounded = 0.0 <= theta_val && theta_val <= 1.0;

  auto log1m_theta = log1p(-theta_val);

  auto n_times_log_theta = elt_multiply(n, log(theta_val));
  auto n_is_zero = n == 0;
  auto N_is_zero = N == 0;
  auto n_is_N = n == N;
  auto N_minus_n = N - n;
  auto logp1 = select(
      N_is_zero, 0.0,
      select(n_is_zero, elt_multiply(N, log1m_theta),
             select(n_is_N, n_times_log_theta,
                    n_times_log_theta + elt_multiply(N_minus_n, log1m_theta))));
  auto logp_expr = colwise_sum(static_select<include_summand<propto>::value>(
      logp1 + binomial_coefficient_log(N, n), logp1));

  auto sum_n_expr = colwise_sum(n);
  auto sum_N_expr = colwise_sum(N);

  auto n_div_theta = elt_divide(n, theta_val);
  auto one_m_theta = 1.0 - theta_val;
  auto deriv_theta = select(
      N_is_zero, 0.0,
      select(n_is_zero, -elt_divide(N, one_m_theta),
             select(n_is_N, n_div_theta,
                    n_div_theta - elt_divide(N_minus_n, one_m_theta))));

  matrix_cl<int> sum_n_cl;
  matrix_cl<int> sum_N_cl;
  matrix_cl<double> logp_cl;
  matrix_cl<double> deriv_cl;

  constexpr bool need_sums
      = !is_constant_all<T_prob_cl>::value && is_stan_scalar<T_prob_cl>::value;
  constexpr bool need_deriv
      = !is_constant_all<T_prob_cl>::value && !is_stan_scalar<T_prob_cl>::value;

  results(check_n_bounded, check_N_nonnegative, check_theta_bounded, logp_cl,
          sum_n_cl, sum_N_cl, deriv_cl)
      = expressions(n_bounded, N_nonnegative, theta_bounded, logp_expr,
                    calc_if<need_sums>(sum_n_expr),
                    calc_if<need_sums>(sum_N_expr),
                    calc_if<need_deriv>(deriv_theta));

  T_partials_return logp = sum(from_matrix_cl(logp_cl));
  auto ops_partials = make_partials_propagator(theta_col);

  if (!is_constant_all<T_prob_cl>::value) {
    if (need_sums) {
      int sum_n = sum(from_matrix_cl(sum_n_cl));
      int sum_N = sum(from_matrix_cl(sum_N_cl));
      double theta_dbl = forward_as<double>(theta_val);
      double& partial = forward_as<internal::broadcast_array<double>>(
          partials<0>(ops_partials))[0];
      if (sum_N != 0) {
        if (sum_n == 0) {
          partial = -sum_N / (1.0 - theta_dbl);
        } else if (sum_n == sum_N) {
          partial = sum_n / theta_dbl;
        } else {
          partial = sum_n / theta_dbl - (sum_N - sum_n) / (1.0 - theta_dbl);
        }
      }
    } else {
      partials<0>(ops_partials) = std::move(deriv_cl);
    }
  }

  return ops_partials.build(logp);
}

}  // namespace math
}  // namespace stan
#endif
#endif
