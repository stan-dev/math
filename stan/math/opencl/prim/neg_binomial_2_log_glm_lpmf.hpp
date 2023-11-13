#ifndef STAN_MATH_OPENCL_PRIM_NEG_BINOMIAL_2_LOG_GLM_LPMF_HPP
#define STAN_MATH_OPENCL_PRIM_NEG_BINOMIAL_2_LOG_GLM_LPMF_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/prim/size.hpp>
#include <stan/math/opencl/rev/operands_and_partials.hpp>
#include <stan/math/opencl/copy.hpp>
#include <stan/math/opencl/prim/multiply.hpp>
#include <stan/math/opencl/plain_type.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/kernels/neg_binomial_2_log_glm_lpmf.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/eval.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/multiply_log.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <vector>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup opencl
 * Returns the log PMF of the Generalized Linear Model (GLM)
 * with Negative-Binomial-2 distribution and log link function.
 * The idea is that neg_binomial_2_log_glm_lpmf(y, x, alpha, beta, phi) should
 * compute a more efficient version of
 * neg_binomial_2_log_lpmf(y, alpha + x * beta, phi) by using analytically
 * simplified gradients.
 * If containers are supplied, returns the log sum of the probabilities.
 * This is an overload of the GLM in
 * prim/prob/neg_binomial_2_log_glm_lpdf.hpp that is implemented in OpenCL.
 * @tparam T_y_cl type of independent variable;
 * this can be a `matrix_cl` vector of intercepts or a single
 * value (wich will be broadcast - used for all instances);
 * @tparam T_x_cl type of the design matrix
 * @tparam T_alpha_cl type of the intercept(s);
 * this can be a vector (of the same length as y) of intercepts or a single
 * value (for models with constant intercept);
 * @tparam T_beta_cl type of the weight vector;
 * this can also be a scalar;
 * @tparam T_phi_cl type of the (positive) precision(s);
 * this can be a vector (of the same length as y, for heteroskedasticity)
 * or a scalar.
 * @param y failures count scalar or vector parameter on OpenCL device. If it
 * is a scalar it will be broadcast - used for all instances.
 * @param x design matrix on OpenCL device. This overload does not support
 * broadcasting of a row vector x!
 * @param alpha intercept (in log odds)
 * @param beta weight vector
 * @param phi (vector of) precision parameter(s)
 * @return log probability or log sum of probabilities
 * @throw std::invalid_argument if container sizes mismatch.
 * @throw std::domain_error if x, beta or alpha is infinite.
 * @throw std::domain_error if phi is infinite or non-positive.
 * @throw std::domain_error if y is negative.
 */
template <bool propto, typename T_y_cl, typename T_x_cl, typename T_alpha_cl,
          typename T_beta_cl, typename T_phi_cl,
          require_all_prim_or_rev_kernel_expression_t<
              T_x_cl, T_y_cl, T_alpha_cl, T_beta_cl, T_phi_cl>* = nullptr>
return_type_t<T_x_cl, T_alpha_cl, T_beta_cl, T_phi_cl>
neg_binomial_2_log_glm_lpmf(const T_y_cl& y, const T_x_cl& x,
                            const T_alpha_cl& alpha, const T_beta_cl& beta,
                            const T_phi_cl& phi) {
  static const char* function = "neg_binomial_2_log_glm_lpmf(OpenCL)";
  using T_partials_return
      = partials_return_t<T_x_cl, T_alpha_cl, T_beta_cl, T_phi_cl>;
  constexpr bool is_y_vector = !is_stan_scalar<T_y_cl>::value;
  constexpr bool is_phi_vector = !is_stan_scalar<T_phi_cl>::value;
  constexpr bool is_alpha_vector = !is_stan_scalar<T_alpha_cl>::value;
  using Eigen::Dynamic;
  using std::isfinite;

  const size_t N = x.rows();
  const size_t M = x.cols();

  if (is_y_vector) {
    check_size_match(function, "Rows of ", "x", N, "rows of ", "y",
                     math::size(y));
  }
  check_size_match(function, "Columns of ", "x", M, "size of ", "beta",
                   math::size(beta));
  if (is_phi_vector) {
    check_size_match(function, "Rows of ", "x", N, "size of ", "phi",
                     math::size(phi));
  }
  if (is_alpha_vector) {
    check_size_match(function, "Rows of ", "x", N, "size of ", "alpha",
                     math::size(alpha));
  }
  if (N == 0) {
    return 0;
  }
  if (!include_summand<propto, T_x_cl, T_alpha_cl, T_beta_cl,
                       T_phi_cl>::value) {
    return 0;
  }

  const auto& y_val = eval(value_of(y));
  const auto& x_val = eval(value_of(x));
  const auto& alpha_val = eval(value_of(alpha));
  const auto& beta_val = eval(value_of(beta));
  const auto& phi_val = eval(value_of(phi));

  // copy any scalars to device, as this is expected by the kernel
  const auto& y_val_cl = to_matrix_cl(y_val);
  const auto& alpha_val_cl = to_matrix_cl(alpha_val);
  const auto& phi_val_cl = to_matrix_cl(phi_val);

  const int local_size
      = opencl_kernels::neg_binomial_2_log_glm.get_option("LOCAL_SIZE_");
  const int wgs = (N + local_size - 1) / local_size;

  const bool need_theta_derivative
      = !is_constant_all<T_x_cl, T_beta_cl, T_alpha_cl>::value;
  matrix_cl<double> theta_derivative_cl(need_theta_derivative ? N : 0, 1);
  const bool need_theta_derivative_sum
      = need_theta_derivative && !is_alpha_vector;
  matrix_cl<double> theta_derivative_sum_cl(wgs, 1);
  const bool need_phi_derivative_sum = !is_alpha_vector;
  const bool need_phi_derivative
      = !is_constant_all<T_phi_cl>::value || need_phi_derivative_sum;
  matrix_cl<double> phi_derivative_cl(
      need_phi_derivative ? (need_phi_derivative_sum ? wgs : N) : 0, 1);
  const bool need_logp1 = include_summand<propto>::value;
  const bool need_logp2
      = include_summand<propto, T_phi_cl>::value && is_phi_vector;
  const bool need_logp3
      = include_summand<propto, T_x_cl, T_alpha_cl, T_beta_cl>::value;
  const bool need_logp4 = include_summand<propto, T_phi_cl>::value
                          && (is_y_vector || is_phi_vector);
  matrix_cl<double> logp_cl(wgs, 1);

  try {
    opencl_kernels::neg_binomial_2_log_glm(
        cl::NDRange(local_size * wgs), cl::NDRange(local_size), logp_cl,
        theta_derivative_cl, theta_derivative_sum_cl, phi_derivative_cl,
        y_val_cl, x_val, alpha_val_cl, beta_val, phi_val_cl, N, M, is_y_vector,
        is_alpha_vector, is_phi_vector, need_theta_derivative,
        need_theta_derivative_sum, need_phi_derivative, need_phi_derivative_sum,
        need_logp1, need_logp2, need_logp3, need_logp4);
  } catch (const cl::Error& e) {
    check_opencl_error(function, e);
  }

  T_partials_return logp = sum(from_matrix_cl(logp_cl));
  if (!std::isfinite(logp)) {
    results(
        check_cl(function, "Vector of dependent variables", y_val,
                 "nonnegative"),
        check_cl(function, "Intercept", alpha_val, "finite"),
        check_cl(function, "Precision parameter", phi_val, "positive finite"))
        = expressions(y_val >= 0, isfinite(alpha_val),
                      isfinite(phi_val) && phi_val > 0);
    check_cl(function, "Design matrix", x_val, "finite") = isfinite(x_val);
    check_cl(function, "Weight vector", beta_val, "finite")
        = isfinite(beta_val);
  } else {
    check_cl(function, "Precision parameter", phi_val, "positive finite")
        = isfinite(phi_val) && phi_val > 0;
  }

  if (include_summand<propto, T_phi_cl>::value && !is_phi_vector) {
    logp += N
            * (multiply_log(forward_as<double>(phi_val),
                            forward_as<double>(phi_val))
               - lgamma(forward_as<double>(phi_val)));
  }
  if (include_summand<propto, T_phi_cl>::value && !is_y_vector
      && !is_phi_vector) {
    logp += forward_as<double>(lgamma(y_val + phi_val)) * N;
  }

  auto ops_partials = make_partials_propagator(x, alpha, beta, phi);
  // Compute the necessary derivatives.
  if (!is_constant<T_x_cl>::value) {
    partials<0>(ops_partials)
        = transpose(beta_val * transpose(theta_derivative_cl));
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
  if (!is_constant_all<T_alpha_cl>::value) {
    if (is_alpha_vector) {
      partials<1>(ops_partials) = std::move(theta_derivative_cl);
    } else {
      forward_as<internal::broadcast_array<double>>(
          partials<1>(ops_partials))[0]
          = sum(from_matrix_cl(theta_derivative_sum_cl));
    }
  }
  if (!is_constant_all<T_phi_cl>::value) {
    if (is_phi_vector) {
      partials<3>(ops_partials) = std::move(phi_derivative_cl);
    } else {
      forward_as<internal::broadcast_array<double>>(
          partials<3>(ops_partials))[0]
          = sum(from_matrix_cl(phi_derivative_cl));
    }
  }
  return ops_partials.build(logp);
}

}  // namespace math
}  // namespace stan

#endif
#endif
