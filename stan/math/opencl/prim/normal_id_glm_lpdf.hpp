#ifndef STAN_MATH_OPENCL_PRIM_NORMAL_ID_GLM_LPDF_HPP
#define STAN_MATH_OPENCL_PRIM_NORMAL_ID_GLM_LPDF_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/elt_divide.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <stan/math/prim/prob/normal_id_glm_lpdf.hpp>
#include <stan/math/opencl/copy.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/multiply.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/rev/operands_and_partials.hpp>
#include <stan/math/opencl/prim/size.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup opencl
 * Returns the log PDF of the Generalized Linear Model (GLM)
 * with Normal distribution and id link function.
 * If containers are supplied, returns the log sum of the probabilities.
 * This is an overload of the GLM in prim/prob/normal_id_glm_lpdf.hpp
 * that is implemented in OpenCL.
 * @tparam T_alpha type of the intercept(s);
 * this can be a (optionally `var_value` containing) `matrix_cl` column vector
 * (of the same length as y) of intercepts or a scalar (for models with
 * constant intercept)
 * @tparam T_beta type of the weight vector;
 * (optionally `var_value` containing) `matrix_cl` column vector
 * @tparam T_scale type of the (positive) scale(s);
 * (optionally `var_value` containing) `matrix_cl` column vector (of the same
 * length as y, for heteroskedasticity) or a scalar.
 * @param y scalar or vector parameter on OpenCL device. If it is a scalar it
 * will be broadcast - used for all instances.
 * @param x design matrix on OpenCL device. This overload does not support
 * broadcasting of a row vector x!
 * @param alpha intercept (in log odds)
 * @param beta weight vector
 * @param sigma (Sequence of) scale parameters for the normal
 * distribution.
 * @return log probability or log sum of probabilities
 * @throw std::domain_error if x, beta or alpha is infinite.
 * @throw std::domain_error if the scale is not positive.
 * @throw std::invalid_argument if container sizes mismatch.
 */
template <bool propto, typename T_x, typename T_y, typename T_alpha,
          typename T_beta, typename T_scale,
          require_all_prim_or_rev_kernel_expression_t<T_x, T_y, T_alpha, T_beta,
                                                      T_scale>* = nullptr>
return_type_t<T_alpha, T_beta, T_scale> normal_id_glm_lpdf(
    const T_y& y, const T_x& x, const T_alpha& alpha, const T_beta& beta,
    const T_scale& sigma) {
  using T_partials_return = partials_return_t<T_alpha, T_beta, T_scale>;
  constexpr bool is_y_vector = is_matrix_cl<decltype(value_of(y))>::value;
  constexpr bool is_sigma_vector
      = is_matrix_cl<decltype(value_of(sigma))>::value;
  constexpr bool is_alpha_vector
      = is_matrix_cl<decltype(value_of(alpha))>::value;
  using T_scale_val = typename std::conditional_t<
      is_sigma_vector, Eigen::Array<partials_return_t<T_scale>, -1, 1>,
      partials_return_t<T_scale>>;
  using Eigen::Array;
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::VectorXd;
  using std::isfinite;
  static const char* function = "normal_id_glm_lpdf(OpenCL)";

  const size_t N = x.rows();
  const size_t M = x.cols();

  if (size(y) != 1) {
    check_size_match(function, "Rows of ", "x_cl", N, "rows of ", "y_cl",
                     size(y));
  }
  check_consistent_size(function, "Weight vector", beta, M);
  if (is_sigma_vector) {
    check_size_match(function, "Rows of ", "x_cl", N, "size of ", "sigma",
                     size(sigma));
  }
  if (is_alpha_vector) {
    check_size_match(function, "Rows of ", "x_cl", N, "size of ", "alpha",
                     size(alpha));
  }
  if (!include_summand<propto, T_alpha, T_beta, T_scale>::value) {
    return 0;
  }
  if (N == 0) {
    return 0;
  }

  const auto& y_val = value_of(y);
  const auto& x_val = value_of(x);
  const auto& alpha_val = value_of(alpha);
  const auto& beta_val = value_of(beta);
  const auto& sigma_val = value_of(sigma);

  auto inv_sigma_expr = elt_divide(1., sigma_val);
  auto y_scaled_expr = elt_multiply(
      (y_val - matrix_vector_multiply(x_val, beta_val) - alpha_val),
      inv_sigma_expr);
  auto mu_derivative_expr = elt_multiply(y_scaled_expr, inv_sigma_expr);
  auto mu_derivative_sum_expr = colwise_sum(mu_derivative_expr);
  auto y_scaled_sq_expr = elt_multiply(y_scaled_expr, y_scaled_expr);
  auto y_scaled_sq_sum_expr = colwise_sum(y_scaled_sq_expr);
  auto sigma_derivative_expr
      = elt_multiply((y_scaled_sq_expr - 1), inv_sigma_expr);
  auto log_sigma_sum_expr = colwise_sum(log(sigma_val));

  const int wgs = y_scaled_sq_sum_expr.rows();

  constexpr bool need_mu_derivative
      = !is_constant<T_beta>::value
        || (!is_constant<T_alpha>::value && is_alpha_vector)
        || (!is_constant<T_y>::value && is_y_vector);
  matrix_cl<double> mu_derivative_cl(need_mu_derivative ? N : 0, 1);
  constexpr bool need_mu_derivative_sum
      = (!is_constant<T_alpha>::value && !is_alpha_vector)
        || (!is_constant<T_y>::value && !is_y_vector);
  matrix_cl<double> mu_derivative_sum_cl(need_mu_derivative_sum ? wgs : 0, 1);
  matrix_cl<double> y_scaled_sq_sum_cl(wgs, 1);
  constexpr bool need_sigma_derivative = !is_constant_all<T_scale>::value;
  matrix_cl<double> sigma_derivative_cl(need_sigma_derivative ? N : 0, 1);
  constexpr bool need_log_sigma_sum
      = include_summand<propto, T_scale>::value && is_sigma_vector;
  matrix_cl<double> log_sigma_sum_cl(need_log_sigma_sum ? wgs : 0, 1);

  results(mu_derivative_cl, mu_derivative_sum_cl, y_scaled_sq_sum_cl,
          sigma_derivative_cl, log_sigma_sum_cl)
      = expressions(calc_if<need_mu_derivative>(mu_derivative_expr),
                    calc_if<need_mu_derivative_sum>(mu_derivative_sum_expr),
                    y_scaled_sq_sum_expr,
                    calc_if<need_sigma_derivative>(sigma_derivative_expr),
                    calc_if<need_log_sigma_sum>(log_sigma_sum_expr));

  double y_scaled_sq_sum = sum(from_matrix_cl(y_scaled_sq_sum_cl));
  operands_and_partials<T_y, T_x, T_alpha, T_beta, T_scale> ops_partials(
      y, x, alpha, beta, sigma);
  double mu_derivative_sum;
  if (need_mu_derivative_sum) {
    mu_derivative_sum = sum(from_matrix_cl(mu_derivative_sum_cl));
  }
  if (!is_constant_all<T_y>::value) {
    if (is_y_vector) {
      ops_partials.edge1_.partials_ = -mu_derivative_cl;
    } else {
      forward_as<internal::broadcast_array<double>>(
          ops_partials.edge1_.partials_)[0]
          = -mu_derivative_sum;
    }
  }
  if (!is_constant_all<T_x>::value) {
    ops_partials.edge2_.partials_
        = transpose(beta_val * transpose(mu_derivative_cl));
  }
  if (!is_constant_all<T_alpha>::value) {
    if (is_alpha_vector) {
      ops_partials.edge3_.partials_ = mu_derivative_cl;
    } else {
      forward_as<internal::broadcast_array<double>>(
          ops_partials.edge3_.partials_)[0]
          = mu_derivative_sum;
    }
  }
  if (!is_constant_all<T_beta>::value) {
    // transposition of a vector can be done without copying
    const matrix_cl<double> mu_derivative_transpose_cl(
        mu_derivative_cl.buffer(), 1, mu_derivative_cl.rows());
    matrix_cl<double>& edge4_partials
        = forward_as<matrix_cl<double>&>(ops_partials.edge4_.partials_);
    matrix_cl<double> edge4_partials_transpose_cl
        = mu_derivative_transpose_cl * x_val;
    edge4_partials = matrix_cl<double>(edge4_partials_transpose_cl.buffer(),
                                       edge4_partials_transpose_cl.cols(), 1);
    if (beta_val.rows() != 0) {
      edge4_partials.add_write_event(
          edge4_partials_transpose_cl.write_events().back());
    }
  }
  if (!is_constant_all<T_scale>::value) {
    ops_partials.edge5_.partials_ = sigma_derivative_cl;
  }

  if (!std::isfinite(y_scaled_sq_sum)) {
    results(
        check_cl(function, "Vector of dependent variables", y_val, "finite"),
        check_cl(function, "Intercept", alpha_val, "finite"),
        check_cl(function, "Scale vector", sigma_val, "positive finite"))
        = expressions(isfinite(y_val), isfinite(alpha_val),
                      isfinite(sigma_val) && sigma_val > 0);
    check_cl(function, "Weight vector", x_val, "finite") = isfinite(x_val);
    check_cl(function, "Weight vector", beta_val, "finite")
        = isfinite(beta_val);
  }

  // Compute log probability.
  T_partials_return logp(0.0);
  if (include_summand<propto>::value) {
    logp += NEG_LOG_SQRT_TWO_PI * N;
  }
  if (include_summand<propto, T_scale>::value) {
    if (is_sigma_vector) {
      logp -= sum(from_matrix_cl(log_sigma_sum_cl));
    } else {
      logp -= N * log(forward_as<double>(sigma_val));
    }
  }
  logp -= 0.5 * y_scaled_sq_sum;

  return ops_partials.build(logp);
}

}  // namespace math
}  // namespace stan

#endif
#endif
