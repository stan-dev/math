#ifndef STAN_MATH_OPENCL_PRIM_POISSON_LOG_GLM_LPMF_HPP
#define STAN_MATH_OPENCL_PRIM_POISSON_LOG_GLM_LPMF_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/prim/size.hpp>
#include <stan/math/opencl/rev/operands_and_partials.hpp>
#include <stan/math/opencl/copy.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/prim/multiply.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup opencl
 * Returns the log PMF of the Generalized Linear Model (GLM)
 * with Poisson distribution and log link function.
 * This is an overload of the GLM in prim/prob/poisson_log_glm_lpmf.hpp
 * that is implemented in OpenCL.
 * @tparam T_y_cl type of independent variable;
 * this can be a `matrix_cl` vector of intercepts or a single
 * value (wich will be broadcast - used for all instances);
 * @tparam T_x_cl type of the design matrix
 * @tparam T_alpha_cl type of the intercept(s);
 * this can be a `matrix_cl` vector (of the same length as y) of intercepts or a
 * single value (for models with constant intercept);
 * @tparam T_beta_cl type of the weight vector;
 * this can also be a single value;
 * @param y positive integer scalar or vector parameter on OpenCL device. If
 * it is a scalar it will be broadcast - used for all instances.
 * @param x design matrix on OpenCL device. This overload does not support
 * broadcasting of a row vector x!
 * @param alpha intercept (in log odds)
 * @param beta weight vector
 * @return log probability or log sum of probabilities
 * @throw std::domain_error if x, beta or alpha is infinite.
 * @throw std::domain_error if y is negative.
 * @throw std::invalid_argument if container sizes mismatch.
 */
template <bool propto, typename T_y_cl, typename T_x_cl, typename T_alpha_cl,
          typename T_beta_cl,
          require_all_prim_or_rev_kernel_expression_t<
              T_y_cl, T_x_cl, T_alpha_cl, T_beta_cl>* = nullptr>
return_type_t<T_x_cl, T_alpha_cl, T_beta_cl> poisson_log_glm_lpmf(
    const T_y_cl& y, const T_x_cl& x, const T_alpha_cl& alpha,
    const T_beta_cl& beta) {
  static const char* function = "poisson_log_glm_lpmf(OpenCL)";
  using T_partials_return = partials_return_t<T_x_cl, T_alpha_cl, T_beta_cl>;
  constexpr bool is_y_vector = !is_stan_scalar<T_y_cl>::value;
  constexpr bool is_alpha_vector = !is_stan_scalar<T_alpha_cl>::value;
  using Eigen::Dynamic;
  using std::exp;
  using std::isfinite;

  const size_t N = x.rows();
  const size_t M = x.cols();

  if (is_y_vector) {
    check_size_match(function, "Rows of ", "x", N, "rows of ", "y",
                     math::size(y));
  }
  check_size_match(function, "Columns of ", "x_cl", M, "size of ", "beta",
                   math::size(beta));
  if (is_alpha_vector) {
    check_size_match(function, "Rows of ", "x", N, "size of ", "alpha",
                     math::size(alpha));
  }
  if (N == 0) {
    return 0;
  }

  if (!include_summand<propto, T_x_cl, T_alpha_cl, T_beta_cl>::value) {
    return 0;
  }

  const auto& y_val = value_of(y);
  const auto& x_val = value_of(x);
  const auto& alpha_val = value_of(alpha);
  const auto& beta_val = value_of(beta);

  T_partials_return logp(0);

  const bool need_logp = include_summand<propto>::value;

  auto theta_expr = matrix_vector_multiply(x_val, beta_val) + alpha_val;
  auto exp_theta_expr = exp(theta_expr);
  auto theta_derivative_expr = select(y_val < 0 || !isfinite(theta_expr),
                                      NOT_A_NUMBER, y_val - exp_theta_expr);
  auto logp_expr
      = colwise_sum(select(need_logp, -lgamma(y_val + 1.0), 0.0)
                    + elt_multiply(y_val, theta_expr) - exp_theta_expr);

  const int wgs = logp_expr.rows();

  matrix_cl<double> theta_derivative_cl(N, 1);
  matrix_cl<double> theta_derivative_sum_cl(wgs, 1);
  matrix_cl<double> logp_cl(wgs, 1);

  results(theta_derivative_cl, theta_derivative_sum_cl, logp_cl) = expressions(
      theta_derivative_expr, colwise_sum(theta_derivative_expr), logp_expr);

  double theta_derivative_sum = sum(from_matrix_cl(theta_derivative_sum_cl));
  logp += sum(from_matrix_cl(logp_cl));
  if (!std::isfinite(theta_derivative_sum)) {
    results(check_cl(function, "Vector of dependent variables", y_val,
                     "nonnegative"),
            check_cl(function, "Intercept", alpha_val, "finite"))
        = expressions(0 <= y_val, isfinite(alpha_val));
    check_cl(function, "Weight vector", beta_val, "finite")
        = isfinite(beta_val);
    check_cl(function, "Matrix of independent variables", x_val, "finite")
        = isfinite(x_val);
  }

  auto ops_partials = make_partials_propagator(x, alpha, beta);
  // Compute the necessary derivatives.
  if (!is_constant_all<T_x_cl>::value) {
    partials<0>(ops_partials)
        = transpose(beta_val * transpose(theta_derivative_cl));
  }
  if (!is_constant_all<T_alpha_cl>::value) {
    if (is_alpha_vector) {
      partials<1>(ops_partials) = theta_derivative_cl;
    } else {
      forward_as<internal::broadcast_array<double>>(
          partials<1>(ops_partials))[0]
          = theta_derivative_sum;
    }
  }
  if (!is_constant_all<T_beta_cl>::value) {
    // transposition of a vector can be done without copying
    const matrix_cl<double> theta_derivative_transpose_cl(
        theta_derivative_cl.buffer(), 1, theta_derivative_cl.rows());
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
