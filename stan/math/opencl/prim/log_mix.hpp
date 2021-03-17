#ifndef STAN_MATH_OPENCL_PRIM_LOG_MIX_HPP
#define STAN_MATH_OPENCL_PRIM_LOG_MIX_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/rev/operands_and_partials.hpp>

namespace stan {
namespace math {

/**
 * Return the log mixture density with specified mixing proportions
 * and log densities.
 *
 * \f[
 * \frac{\partial }{\partial p_x}
 * \log\left(\exp^{\log\left(p_1\right)+d_1}+\cdot\cdot\cdot+
 * \exp^{\log\left(p_n\right)+d_n}\right)
 * =\frac{e^{d_x}}{e^{d_1}p_1+\cdot\cdot\cdot+e^{d_m}p_m}
 * \f]
 *
 * \f[
 * \frac{\partial }{\partial d_x}
 * \log\left(\exp^{\log\left(p_1\right)+d_1}+\cdot\cdot\cdot+
 * \exp^{\log\left(p_n\right)+d_n}\right)
 * =\frac{e^{d_x}p_x}{e^{d_1}p_1+\cdot\cdot\cdot+e^{d_m}p_m}
 * \f]
 *
 * @tparam T_theta Type of theta.
 * @tparam T_lam Type of lambda.
 * @param theta std/row/col vector of mixing proportions in [0, 1].
 * @param lambda std/row/col vector of log densities.
 * @return log mixture of densities in specified proportion
 */
template <typename T_theta_cl, typename T_lambda_cl,
          require_all_prim_or_rev_kernel_expression_t<T_theta_cl,
                                                      T_lambda_cl>* = nullptr>
inline auto log_mix(const T_theta_cl& theta, const T_lambda_cl& lambda) {
  static const char* function = "log_mix(OpenCL)";
  using T_return = return_type_t<T_theta_cl, T_lambda_cl>;
  using std::isfinite;
  using std::isnan;

  check_consistent_sizes(function, "theta", theta, "lambda", lambda);
  const size_t N = max_size(theta, lambda);
  if (N == 0) {
    return T_return(0.0);
  }

  const auto& theta_col = as_column_vector_or_scalar(theta);

  const auto& theta_val = value_of(theta_col);
  const auto& lambda_val = value_of(lambda);

  auto check_lambda_not_nan
      = check_cl(function, "lambda", lambda_val, "not NaN");
  auto lambda_not_nan = !isnan(lambda_val);
  auto check_theta_bounded
      = check_cl(function, "theta", theta_val, "in the interval[0, 1]");
  auto theta_bounded = 0.0 <= theta_val && theta_val <= 1.0;

  auto theta_bc = rowwise_broadcast(theta_val);
  auto lambda_p_log_theta_expr = lambda_val + log(theta_bc);
  matrix_cl<double> lambda_p_log_theta;
  matrix_cl<double> lambda_p_log_theta_colwise_max;
  if (theta.cols() == lambda.cols()) {
    results(check_lambda_not_nan, check_theta_bounded, lambda_p_log_theta,
            lambda_p_log_theta_colwise_max)
        = expressions(lambda_not_nan, theta_bounded, lambda_p_log_theta_expr,
                      colwise_max(lambda_p_log_theta_expr));
  } else {
    results(check_lambda_not_nan, lambda_p_log_theta,
            lambda_p_log_theta_colwise_max)
        = expressions(lambda_not_nan, lambda_p_log_theta_expr,
                      colwise_max(lambda_p_log_theta_expr));
    check_theta_bounded = theta_bounded;
  }
  while (lambda_p_log_theta_colwise_max.rows() > 1) {
    lambda_p_log_theta_colwise_max
        = colwise_max(lambda_p_log_theta_colwise_max).eval();
  }
  matrix_cl<double> sum_exp = colwise_sum(exp(
      lambda_p_log_theta - colwise_broadcast(lambda_p_log_theta_colwise_max)));
  while (sum_exp.rows() > 1) {
    sum_exp = colwise_sum(sum_exp).eval();
  }

  auto logp_vec_expr = transpose(lambda_p_log_theta_colwise_max + log(sum_exp));
  matrix_cl<double> logp_vec;
  matrix_cl<double> logp_sum;
  results(logp_vec, logp_sum) = expressions(
      calc_if<!is_constant_all<T_theta_cl, T_lambda_cl>::value>(logp_vec_expr),
      colwise_sum(logp_vec_expr));

  operands_and_partials<decltype(theta_col), decltype(lambda)> ops_partials(
      theta_col, lambda);
  if (!is_constant_all<T_theta_cl, T_lambda_cl>::value) {
    auto derivs_expr = exp(lambda_val - colwise_broadcast(transpose(logp_vec)));
    if (!is_constant<T_lambda_cl>::value) {
      auto lambda_deriv_expr = elt_multiply(derivs_expr, theta_bc);
      matrix_cl<double> derivs;
      matrix_cl<double> lambda_deriv;
      results(derivs, lambda_deriv)
          = expressions(calc_if<!is_constant<T_theta_cl>::value>(derivs_expr),
                        lambda_deriv_expr);

      ops_partials.edge2_.partials_ = std::move(lambda_deriv);
      if (!is_constant<T_theta_cl>::value) {
        ops_partials.edge1_.partials_ = rowwise_sum(derivs);
      }
    } else if (!is_constant<T_theta_cl>::value) {
      ops_partials.edge1_.partials_ = rowwise_sum(derivs_expr);
    }
  }
  return ops_partials.build(sum(from_matrix_cl(logp_sum)));
}
}  // namespace math
}  // namespace stan
#endif
#endif
