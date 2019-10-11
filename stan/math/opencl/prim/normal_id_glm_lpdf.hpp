#ifndef STAN_MATH_OPENCL_PRIM_NORMAL_ID_GLM_LPDF_HPP
#define STAN_MATH_OPENCL_PRIM_NORMAL_ID_GLM_LPDF_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/mat/fun/value_of_rec.hpp>
#include <stan/math/prim/scal/fun/size_zero.hpp>
#include <stan/math/prim/scal/fun/sum.hpp>
#include <stan/math/prim/arr/fun/value_of_rec.hpp>
#include <stan/math/prim/mat/prob/normal_id_glm_lpdf.hpp>

#include <stan/math/opencl/kernels/normal_id_glm_lpdf.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/prim/multiply.hpp>

#include <cmath>

namespace stan {
namespace math {

/**
 * Returns the log PDF of the Generalized Linear Model (GLM)
 * with Normal distribution and id link function.
 * If containers are supplied, returns the log sum of the probabilities.
 * This is an overload of the GLM in prim/mar/prob/normal_id_glm_lpdf.hpp
 * that is implemented in OpenCL.
 * @tparam T_alpha type of the intercept(s);
 * this can be a vector (of the same length as y) of intercepts or a single
 * value (for models with constant intercept);
 * @tparam T_beta type of the weight vector;
 * this can also be a single value;
 * @tparam T_scale type of the (positive) scale(s);
 * this can be a vector (of the same length as y, for heteroskedasticity)
 * or a scalar.
 * @param y_cl vector parameter on OpenCL device
 * @param x_cl design matrix on OpenCL device
 * @param alpha intercept (in log odds)
 * @param beta weight vector
 * @param sigma (Sequence of) scale parameters for the normal
 * distribution.
 * @return log probability or log sum of probabilities
 * @throw std::domain_error if x, beta or alpha is infinite.
 * @throw std::domain_error if the scale is not positive.
 * @throw std::invalid_argument if container sizes mismatch.
 */
template <bool propto, typename T_alpha, typename T_beta, typename T_scale>
return_type_t<T_alpha, T_beta, T_scale> normal_id_glm_lpdf(
    const matrix_cl<double> &y_cl, const matrix_cl<double> &x_cl,
    const T_alpha &alpha, const T_beta &beta, const T_scale &sigma) {
  static const char *function = "normal_id_glm_lpdf(OpenCL)";
  using T_partials_return = partials_return_t<T_alpha, T_beta, T_scale>;
  using T_scale_val = typename std::conditional_t<
      is_vector<T_scale>::value,
      Eigen::Array<partials_return_t<T_scale>, -1, 1>,
      partials_return_t<T_scale>>;

  using Eigen::Array;
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::VectorXd;

  const size_t N = x_cl.rows();
  const size_t M = x_cl.cols();

  check_positive_finite(function, "Scale vector", sigma);
  check_size_match(function, "Rows of ", "x_cl", N, "rows of ", "y_cl",
                   y_cl.rows());
  check_consistent_size(function, "Weight vector", beta, M);
  if (is_vector<T_scale>::value) {
    check_size_match(function, "Rows of ", "y_cl", N, "size of ", "sigma",
                     length(sigma));
  }
  if (is_vector<T_alpha>::value) {
    check_size_match(function, "Rows of ", "y_cl", N, "size of ", "alpha",
                     length(alpha));
  }

  if (!include_summand<propto, T_alpha, T_beta, T_scale>::value) {
    return 0;
  }

  if (N == 0) {
    return 0;
  }

  const auto &beta_val = value_of_rec(beta);
  const auto &alpha_val = value_of_rec(alpha);
  const auto &sigma_val = value_of_rec(sigma);

  const auto &beta_val_vec = as_column_vector_or_scalar(beta_val);
  const auto &alpha_val_vec = as_column_vector_or_scalar(alpha_val);
  const auto &sigma_val_vec = as_column_vector_or_scalar(sigma_val);

  T_scale_val inv_sigma = 1 / as_array_or_scalar(sigma_val_vec);
  Matrix<T_partials_return, Dynamic, 1> y_minus_mu_over_sigma_mat(N);
  auto y_scaled = y_minus_mu_over_sigma_mat.array();

  const int local_size
      = opencl_kernels::normal_id_glm.get_option("LOCAL_SIZE_");
  const int wgs = (N + local_size - 1) / local_size;

  matrix_cl<double> beta_cl(beta_val_vec);
  matrix_cl<double> alpha_cl(alpha_val_vec);
  matrix_cl<double> sigma_cl(sigma_val_vec);

  const bool need_mu_derivative = !is_constant_all<T_beta, T_alpha>::value;
  matrix_cl<double> mu_derivative_cl(need_mu_derivative ? N : 0, 1);
  const bool need_mu_derivative_sum
      = !is_constant_all<T_alpha>::value && !is_vector<T_alpha>::value;
  matrix_cl<double> mu_derivative_sum_cl(need_mu_derivative_sum ? wgs : 0, 1);
  matrix_cl<double> y_minus_mu_over_sigma_squared_sum_cl(wgs, 1);
  const bool need_sigma_derivative
      = !is_constant_all<T_scale>::value && is_vector<T_scale>::value;
  matrix_cl<double> sigma_derivative_cl(need_sigma_derivative ? N : 0, 1);
  const bool need_log_sigma_sum
      = include_summand<propto, T_scale>::value && is_vector<T_scale>::value;
  matrix_cl<double> log_sigma_sum_cl(need_log_sigma_sum ? wgs : 0, 1);

  try {
    opencl_kernels::normal_id_glm(
        cl::NDRange(local_size * wgs), cl::NDRange(local_size),
        mu_derivative_cl, mu_derivative_sum_cl,
        y_minus_mu_over_sigma_squared_sum_cl, sigma_derivative_cl,
        log_sigma_sum_cl, y_cl, x_cl, alpha_cl, beta_cl, sigma_cl, N, M,
        length(alpha) != 1, length(sigma) != 1, need_mu_derivative,
        need_mu_derivative_sum, need_sigma_derivative, need_log_sigma_sum);
  } catch (const cl::Error &e) {
    check_opencl_error(function, e);
  }
  double y_scaled_sq_sum
      = sum(from_matrix_cl(y_minus_mu_over_sigma_squared_sum_cl));

  operands_and_partials<T_alpha, T_beta, T_scale> ops_partials(alpha, beta,
                                                               sigma);

  if (!is_constant_all<T_alpha>::value && is_vector<T_alpha>::value) {
    ops_partials.edge1_.partials_
        = from_matrix_cl<Dynamic, 1>(mu_derivative_cl);
  }
  if (need_mu_derivative_sum) {
    ops_partials.edge1_.partials_[0]
        = sum(from_matrix_cl(mu_derivative_sum_cl));
  }
  if (!is_constant_all<T_beta>::value) {
    const matrix_cl<double> mu_derivative_transpose_cl(
        mu_derivative_cl.buffer(), 1,
        mu_derivative_cl.rows());  // transposition of a vector can be done
    // without copying
    ops_partials.edge2_.partials_
        = from_matrix_cl<1, Dynamic>(mu_derivative_transpose_cl * x_cl);
  }
  if (!is_constant_all<T_scale>::value) {
    if (is_vector<T_scale>::value) {
      ops_partials.edge3_.partials_
          = from_matrix_cl<Dynamic, 1>(sigma_derivative_cl);
    } else {
      ops_partials.edge3_.partials_[0]
          = (y_scaled_sq_sum - N) * as_scalar(inv_sigma);
    }
  }

  if (!std::isfinite(y_scaled_sq_sum)) {
    check_finite(function, "Vector of dependent variables",
                 from_matrix_cl(y_cl));
    check_finite(function, "Weight vector", beta);
    check_finite(function, "Intercept", alpha);
    // if all other checks passed, next will only fail if x is not finite
    check_finite(function, "Matrix of independent variables", y_scaled_sq_sum);
  }

  // Compute log probability.
  T_partials_return logp(0.0);
  if (include_summand<propto>::value) {
    logp += NEG_LOG_SQRT_TWO_PI * N;
  }
  if (include_summand<propto, T_scale>::value) {
    if (is_vector<T_scale>::value) {
      logp -= sum(from_matrix_cl(log_sigma_sum_cl));
    } else {
      logp -= N * log(as_scalar(sigma_val));
    }
  }
  if (include_summand<propto, T_alpha, T_beta, T_scale>::value) {
    logp -= 0.5 * y_scaled_sq_sum;
  }
  return ops_partials.build(logp);
}

template <typename T_alpha, typename T_beta, typename T_scale>
inline return_type_t<T_alpha, T_beta, T_scale> normal_id_glm_lpdf(
    const matrix_cl<double> &y, const matrix_cl<double> &x,
    const T_alpha &alpha, const T_beta &beta, const T_scale &sigma) {
  return normal_id_glm_lpdf<false>(y, x, alpha, beta, sigma);
}

}  // namespace math
}  // namespace stan

#endif
#endif
