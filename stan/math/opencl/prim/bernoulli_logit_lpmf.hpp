#ifndef STAN_MATH_OPENCL_PRIM_BERNOULLI_LOGIT_LPMF_HPP
#define STAN_MATH_OPENCL_PRIM_BERNOULLI_LOGIT_LPMF_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Returns the log PMF of the logit-parametrized Bernoulli distribution. If
 * containers are supplied, returns the log sum of the probabilities.
 *
 * @tparam T_n_cl type of integer parameter
 * @tparam T_prob_cl type of chance of success parameter
 * @param n integer parameter
 * @param theta logit-transformed chance of success parameter
 * @return log probability or log sum of probabilities
 * @throw std::domain_error if theta is infinite.
 * @throw std::invalid_argument if container sizes mismatch.
 */
template <
    bool propto, typename T_n_cl, typename T_prob_cl,
    require_all_prim_or_rev_kernel_expression_t<T_n_cl, T_prob_cl>* = nullptr,
    require_any_not_stan_scalar_t<T_n_cl, T_prob_cl>* = nullptr>
return_type_t<T_prob_cl> bernoulli_logit_lpmf(const T_n_cl& n,
                                              const T_prob_cl& theta) {
  static const char* function = "bernoulli_logit_lpmf(OpenCL)";
  using T_partials_return = partials_return_t<T_prob_cl>;
  using std::isnan;
  constexpr bool is_n_vector = !is_stan_scalar<T_n_cl>::value;
  constexpr bool is_theta_vector = !is_stan_scalar<T_prob_cl>::value;

  check_consistent_sizes(function, "Random variable", n,
                         "Probability parameter", theta);
  const size_t N = is_n_vector ? size(n) : size(theta);
  if (N == 0) {
    return 0.0;
  }
  if (!include_summand<propto, T_prob_cl>::value) {
    return 0.0;
  }

  const auto& theta_col = as_column_vector_or_scalar(theta);
  const auto& theta_val = value_of(theta_col);

  auto check_n_bounded = check_cl(function, "n", n, "in the interval [0, 1]");
  auto n_bounded_expr = 0 <= n && n <= 1;
  auto check_theta_not_nan
      = check_cl(function, "Logit transformed probability parameter", theta_val,
                 "not NaN");
  auto theta_not_nan_expr = !isnan(theta_val);

  auto signs_expr = 2 * n - 1.0;  // subtracting 1.0 converts int to double
  auto ntheta_expr = elt_multiply(signs_expr, theta_val);
  auto exp_m_ntheta_expr = exp(-ntheta_expr);
  static const double cutoff = 20.0;
  auto condition1_expr = ntheta_expr > cutoff;
  auto condition2_expr = ntheta_expr < -cutoff;
  auto logp_expr = colwise_sum(
      select(condition1_expr, -exp_m_ntheta_expr,
             select(condition2_expr, ntheta_expr, -log1p(exp_m_ntheta_expr))));
  auto deriv_expr = select(
      condition1_expr, -exp_m_ntheta_expr,
      select(condition2_expr, signs_expr,
             elt_multiply(signs_expr, elt_divide(exp_m_ntheta_expr,
                                                 (exp_m_ntheta_expr + 1)))));

  matrix_cl<double> logp_cl;
  matrix_cl<double> deriv_cl;

  results(logp_cl, deriv_cl, check_n_bounded, check_theta_not_nan)
      = expressions(logp_expr,
                    calc_if<!is_constant_all<T_prob_cl>::value>(deriv_expr),
                    n_bounded_expr, theta_not_nan_expr);

  T_partials_return logp = sum(from_matrix_cl(logp_cl));
  operands_and_partials<decltype(theta_col)> ops_partials(theta_col);

  if (!is_constant_all<T_prob_cl>::value) {
    ops_partials.edge1_.partials_ = deriv_cl;
  }

  return ops_partials.build(logp);
}

}  // namespace math
}  // namespace stan
#endif
#endif
