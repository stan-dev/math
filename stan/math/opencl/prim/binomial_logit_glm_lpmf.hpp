#ifndef STAN_MATH_OPENCL_PRIM_BINOMIAL_LOGIT_GLM_LPMF_HPP
#define STAN_MATH_OPENCL_PRIM_BINOMIAL_LOGIT_GLM_LPMF_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/prim/size.hpp>
#include <stan/math/opencl/rev/operands_and_partials.hpp>
#include <stan/math/opencl/copy.hpp>
#include <stan/math/opencl/prim/multiply.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/max.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>

#include <cmath>
#include <cstdint>

namespace stan {
namespace math {

template <bool propto, typename T_n_cl, typename T_N_cl, typename T_x_cl,
          typename T_alpha_cl, typename T_beta_cl,
          require_all_prim_or_rev_kernel_expression_t<
              T_n_cl, T_N_cl, T_x_cl, T_alpha_cl, T_beta_cl>* = nullptr>
return_type_t<T_x_cl, T_alpha_cl, T_beta_cl> binomial_logit_glm_lpmf(
    const T_n_cl& n, const T_N_cl& N, const T_x_cl& x, const T_alpha_cl& alpha,
    const T_beta_cl& beta) {
  static const char* function = "binomial_logit_glm_lpmf(OpenCL)";
  using T_partials_return = partials_return_t<T_x_cl, T_alpha_cl, T_beta_cl>;
  constexpr bool is_y_vector = !is_stan_scalar<T_n_cl>::value;
  constexpr bool is_alpha_vector = !is_stan_scalar<T_alpha_cl>::value;

  const size_t N_instances
      = max(max_size(n, N, alpha), static_cast<int64_t>(x.rows()));
  const size_t N_attributes = x.cols();

  check_consistent_sizes(function, "Successes variable", n,
                         "Population size parameter", N);
  check_consistent_size(function, "Successes variable", n, N_instances);
  check_consistent_size(function, "Population size parameter", N, N_instances);
  check_consistent_size(function, "Weight vector", beta, N_attributes);
  check_consistent_size(function, "Vector of intercepts", alpha, N_instances);

  if (N_instances == 0 || N_attributes == 0) {
    return 0;
  }
  if (!include_summand<propto, T_x_cl, T_alpha_cl, T_beta_cl>::value) {
    return 0;
  }

  auto&& x_val = value_of(x);
  auto&& alpha_val = value_of(alpha);
  auto&& beta_val = value_of(beta);

  auto check_n_bounded
      = check_cl(function, "Successes variable", n, "in the interval [0, N]");
  auto n_bounded = 0 <= n && n <= N;
  auto check_N_nonnegative
      = check_cl(function, "Population size variable", n, "nonnegative");
  auto N_nonnegative = N >= 0;

  auto theta_expr = matrix_vector_multiply(x_val, beta_val) + alpha_val;
  auto log_inv_logit_theta = log_inv_logit(theta_expr);
  auto log1m_inv_logit_theta = log1m_inv_logit(theta_expr);
  auto n_diff = N - n;
  auto logp_expr1 = elt_multiply(n, log_inv_logit_theta)
                    + elt_multiply(n_diff, log1m_inv_logit_theta);
  auto logp_expr
      = static_select<include_summand<propto, T_n_cl, T_N_cl>::value>(
          logp_expr1 + binomial_coefficient_log(N, n), logp_expr1);

  constexpr bool need_theta_deriv
      = !is_constant_all<T_beta_cl, T_x_cl, T_alpha_cl>::value;
  auto theta_deriv_expr = n - elt_multiply(N, exp(log_inv_logit_theta));

  constexpr bool need_theta_deriv_sum = need_theta_deriv && !is_alpha_vector;
  matrix_cl<double> logp_cl;
  matrix_cl<double> theta_deriv_cl;
  matrix_cl<double> theta_deriv_sum_cl;

  results(check_n_bounded, check_N_nonnegative, logp_cl, theta_deriv_cl,
          theta_deriv_sum_cl)
      = expressions(
          n_bounded, N_nonnegative, logp_expr,
          calc_if<need_theta_deriv>(theta_deriv_expr),
          calc_if<need_theta_deriv_sum>(colwise_sum(theta_deriv_expr)));

  T_partials_return logp = sum(from_matrix_cl(logp_cl));
  using std::isfinite;
  if (!isfinite(logp)) {
    check_cl(function, "Intercept", alpha_val, "finite") = isfinite(alpha_val);
    check_cl(function, "Weight vector", beta_val, "finite")
        = isfinite(beta_val);
    check_cl(function, "Matrix of independent variables", x_val, "finite")
        = isfinite(x_val);
  }

  auto ops_partials = make_partials_propagator(x, alpha, beta);
  if (!is_constant_all<T_x_cl>::value) {
    partials<0>(ops_partials) = transpose(beta_val * transpose(theta_deriv_cl));
  }
  if (!is_constant_all<T_alpha_cl>::value) {
    if (is_alpha_vector) {
      partials<1>(ops_partials) = theta_deriv_cl;
    } else {
      forward_as<internal::broadcast_array<double>>(
          partials<1>(ops_partials))[0]
          = sum(from_matrix_cl(theta_deriv_sum_cl));
    }
  }
  if (!is_constant_all<T_beta_cl>::value) {
    // transposition of a vector can be done without copying
    const matrix_cl<double> theta_derivative_transpose_cl(
        theta_deriv_cl.buffer(), 1, theta_deriv_cl.rows());
    matrix_cl<double> edge3_partials_transpose_cl
        = theta_derivative_transpose_cl * x_val;
    partials<2>(ops_partials)
        = matrix_cl<double>(edge3_partials_transpose_cl.buffer(),
                            edge3_partials_transpose_cl.cols(), 1);
    if (beta_val.rows() != 0) {
      edge<2>(ops_partials)
          .partials_.add_write_event(
              edge3_partials_transpose_cl.write_events().back());
    }
  }
  return ops_partials.build(logp);
}

}  // namespace math
}  // namespace stan

#endif
#endif
