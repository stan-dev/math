#ifndef STAN_MATH_OPENCL_PRIM_NEG_BINOMIAL_2_LOG_GLM_LPMF_HPP
#define STAN_MATH_OPENCL_PRIM_NEG_BINOMIAL_2_LOG_GLM_LPMF_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <stan/math/prim/scal/err/check_nonnegative.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/fun/multiply_log.hpp>
#include <stan/math/prim/scal/fun/digamma.hpp>
#include <stan/math/prim/mat/fun/lgamma.hpp>
#include <stan/math/prim/mat/fun/value_of_rec.hpp>
#include <stan/math/prim/arr/fun/value_of_rec.hpp>
#include <stan/math/prim/scal/fun/sum.hpp>

#include <stan/math/opencl/kernels/neg_binomial_2_log_glm_lpmf.hpp>

#include <vector>
#include <cmath>

namespace stan {
namespace math {

/**
 * Returns the log PMF of the Generalized Linear Model (GLM)
 * with Negative-Binomial-2 distribution and log link function.
 * The idea is that neg_binomial_2_log_glm_lpmf(y, x, alpha, beta, phi) should
 * compute a more efficient version of
 * neg_binomial_2_log_lpmf(y, alpha + x * beta, phi) by using analytically
 * simplified gradients.
 * If containers are supplied, returns the log sum of the probabilities.
 * This is an overload of the GLM in
 * prim/mar/prob/neg_binomial_2_log_glm_lpdf.hpp that is implemented in OpenCL.
 * @tparam T_alpha type of the intercept(s);
 * this can be a vector (of the same length as y) of intercepts or a single
 * value (for models with constant intercept);
 * @tparam T_beta type of the weight vector;
 * this can also be a scalar;
 * @tparam T_precision type of the (positive) precision(s);
 * this can be a vector (of the same length as y, for heteroskedasticity)
 * or a scalar.
 * @param y_cl failures count vector parameter on OpenCL device
 * @param x_cl design matrix on OpenCL device
 * @param alpha intercept (in log odds)
 * @param beta weight vector
 * @param phi (vector of) precision parameter(s)
 * @return log probability or log sum of probabilities
 * @throw std::invalid_argument if container sizes mismatch.
 * @throw std::domain_error if x, beta or alpha is infinite.
 * @throw std::domain_error if phi is infinite or non-positive.
 * @throw std::domain_error if y is negative.
 */
template <bool propto, typename T_alpha, typename T_beta, typename T_precision>
return_type_t<T_alpha, T_beta, T_precision> neg_binomial_2_log_glm_lpmf(
    const matrix_cl<int>& y_cl, const matrix_cl<double>& x_cl,
    const T_alpha& alpha, const T_beta& beta, const T_precision& phi) {
  static const char* function = "neg_binomial_2_log_glm_lpmf(OpenCL)";
  using T_partials_return = partials_return_t<T_alpha, T_beta, T_precision>;

  using Eigen::Array;
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::VectorXd;

  T_partials_return logp(0.0);

  const size_t N = x_cl.rows();
  const size_t M = x_cl.cols();

  check_size_match(function, "Rows of ", "x_cl", N, "rows of ", "y_cl",
                   y_cl.rows());
  check_consistent_size(function, "Weight vector", beta, M);
  if (is_vector<T_precision>::value) {
    check_size_match(function, "Rows of ", "y_cl", N, "size of ", "phi",
                     length(phi));
  }
  if (is_vector<T_alpha>::value) {
    check_size_match(function, "Rows of ", "y_cl", N, "size of ", "alpha",
                     length(alpha));
  }
  check_positive_finite(function, "Precision parameter", phi);

  if (N == 0) {
    return 0;
  }

  if (!include_summand<propto, T_alpha, T_beta, T_precision>::value) {
    return 0;
  }

  const auto& beta_val = value_of_rec(beta);
  const auto& alpha_val = value_of_rec(alpha);
  const auto& phi_val = value_of_rec(phi);

  const auto& beta_val_vec = as_column_vector_or_scalar(beta_val);
  const auto& alpha_val_vec = as_column_vector_or_scalar(alpha_val);
  const auto& phi_val_vec = as_column_vector_or_scalar(phi_val);

  const auto& phi_arr = as_array_or_scalar(phi_val_vec);

  const int local_size
      = opencl_kernels::neg_binomial_2_log_glm.get_option("LOCAL_SIZE_");
  const int wgs = (N + local_size - 1) / local_size;

  const matrix_cl<double> beta_cl(beta_val_vec);
  const matrix_cl<double> alpha_cl(alpha_val_vec);
  const matrix_cl<double> phi_cl(phi_val_vec);

  const bool need_theta_derivative = !is_constant_all<T_beta, T_alpha>::value;
  matrix_cl<double> theta_derivative_cl(need_theta_derivative ? N : 0, 1);
  const bool need_theta_derivative_sum
      = need_theta_derivative && !is_vector<T_alpha>::value;
  matrix_cl<double> theta_derivative_sum_cl(wgs, 1);
  const bool need_phi_derivative_sum = !is_vector<T_alpha>::value;
  const bool need_phi_derivative
      = !is_constant_all<T_precision>::value || need_phi_derivative_sum;
  matrix_cl<double> phi_derivative_cl(
      need_phi_derivative ? (need_phi_derivative_sum ? wgs : N) : 0, 1);
  const bool need_logp1 = include_summand<propto>::value;
  const bool need_logp2 = include_summand<propto, T_precision>::value
                          && is_vector<T_precision>::value;
  const bool need_logp3
      = include_summand<propto, T_alpha, T_beta, T_precision>::value;
  const bool need_logp4 = include_summand<propto, T_alpha, T_beta>::value;
  const bool need_logp5 = include_summand<propto, T_precision>::value;
  const bool need_logp
      = need_logp1 || need_logp2 || need_logp3 || need_logp4 || need_logp5;
  matrix_cl<double> logp_cl(need_logp ? wgs : 0, 1);

  try {
    opencl_kernels::neg_binomial_2_log_glm(
        cl::NDRange(local_size * wgs), cl::NDRange(local_size), logp_cl,
        theta_derivative_cl, theta_derivative_sum_cl, phi_derivative_cl, y_cl,
        x_cl, alpha_cl, beta_cl, phi_cl, N, M, length(alpha) != 1,
        length(phi) != 1, need_theta_derivative, need_theta_derivative_sum,
        need_phi_derivative, need_phi_derivative_sum, need_logp1, need_logp2,
        need_logp3, need_logp4, need_logp5);
  } catch (const cl::Error& e) {
    check_opencl_error(function, e);
  }

  double logp_sum = sum(from_matrix_cl(logp_cl));
  if (!std::isfinite(logp_sum)) {
    check_nonnegative(function, "Failures variables", from_matrix_cl(y_cl));
    check_finite(function, "Weight vector", beta);
    check_finite(function, "Intercept", alpha);
    // if all other checks passed, next will only fail if x is not finite
    check_finite(function, "Matrix of independent variables", logp_sum);
  }
  if (need_logp) {
    logp += logp_sum;
  }

  if (include_summand<propto, T_precision>::value
      && !is_vector<T_precision>::value) {
    logp += N
            * (multiply_log(as_scalar(phi_val), as_scalar(phi_val))
               - lgamma(as_scalar(phi_val)));
  }

  operands_and_partials<T_alpha, T_beta, T_precision> ops_partials(alpha, beta,
                                                                   phi);
  // Compute the necessary derivatives.
  if (!is_constant_all<T_alpha>::value) {
    if (is_vector<T_alpha>::value) {
      ops_partials.edge1_.partials_
          = from_matrix_cl<Dynamic, 1>(theta_derivative_cl);
    } else {
      ops_partials.edge1_.partials_[0]
          = sum(from_matrix_cl<Dynamic, 1>(theta_derivative_sum_cl));
    }
  }
  if (!is_constant_all<T_beta>::value) {
    matrix_cl<double> theta_derivative_transpose_cl(
        theta_derivative_cl.buffer(), 1,
        theta_derivative_cl
            .rows());  // transposition of a vector can be done without copying
    ops_partials.edge2_.partials_
        = from_matrix_cl<1, Dynamic>(theta_derivative_transpose_cl * x_cl);
  }
  if (!is_constant_all<T_precision>::value) {
    if (is_vector<T_precision>::value) {
      ops_partials.edge3_.partials_
          = std::move(from_matrix_cl<Dynamic, 1>(phi_derivative_cl));
    } else {
      ops_partials.edge3_.partials_[0]
          = sum(from_matrix_cl<Dynamic, 1>(phi_derivative_cl));
    }
  }
  return ops_partials.build(logp);
}

template <typename T_alpha, typename T_beta, typename T_precision>
inline return_type_t<T_alpha, T_beta, T_precision> neg_binomial_2_log_glm_lpmf(
    const matrix_cl<int>& y, const matrix_cl<double>& x, const T_alpha& alpha,
    const T_beta& beta, const T_precision& phi) {
  return neg_binomial_2_log_glm_lpmf<false>(y, x, alpha, beta, phi);
}

}  // namespace math
}  // namespace stan

#endif
#endif
