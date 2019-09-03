#ifndef STAN_MATH_PRIM_MAT_PROB_NORMAL_ID_GLM_LPDF_HPP
#define STAN_MATH_PRIM_MAT_PROB_NORMAL_ID_GLM_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/mat/fun/value_of_rec.hpp>
#include <stan/math/prim/scal/fun/size_zero.hpp>
#include <stan/math/prim/scal/fun/sum.hpp>
#include <stan/math/prim/arr/fun/value_of_rec.hpp>

#ifdef STAN_OPENCL
#include <stan/math/opencl/kernels/normal_id_glm_lpdf.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/multiply.hpp>
#endif

#include <cmath>

namespace stan {
namespace math {

/**
 * Returns the log PDF of the Generalized Linear Model (GLM)
 * with Normal distribution and id link function.
 * If containers are supplied, returns the log sum of the probabilities.
 * The idea is that normal_id_glm_lpdf(y, x, alpha, beta, sigma) should
 * compute a more efficient version of normal_lpdf(y, alpha + x * beta, sigma)
 * by using analytically simplified gradients.
 * @tparam T_y type of vector of dependent variables (labels);
 * @tparam T_x type of the matrix of independent variables (features); this
 * should be an Eigen::Matrix type whose number of rows should match the
 * length of y and whose number of columns should match the length of beta
 * @tparam T_alpha type of the intercept(s);
 * this can be a vector (of the same length as y) of intercepts or a single
 * value (for models with constant intercept);
 * @tparam T_beta type of the weight vector;
 * this can also be a single value;
 * @tparam T_scale type of the (positive) scale(s);
 * this can be a vector (of the same length as y, for heteroskedasticity)
 * or a scalar.
 * @param y vector parameter
 * @param x design matrix
 * @param alpha intercept (in log odds)
 * @param beta weight vector
 * @param sigma (Sequence of) scale parameters for the normal
 * distribution.
 * @return log probability or log sum of probabilities
 * @throw std::domain_error if x, beta or alpha is infinite.
 * @throw std::domain_error if the scale is not positive.
 * @throw std::invalid_argument if container sizes mismatch.
 */
template <bool propto, typename T_y, typename T_x, typename T_alpha,
          typename T_beta, typename T_scale>
return_type_t<T_y, T_x, T_alpha, T_beta, T_scale> normal_id_glm_lpdf(
    const T_y &y, const T_x &x, const T_alpha &alpha, const T_beta &beta,
    const T_scale &sigma) {
  static const char *function = "normal_id_glm_lpdf";
  typedef partials_return_type_t<T_y, T_x, T_alpha, T_beta, T_scale>
      T_partials_return;
  typedef typename std::conditional_t<
      is_vector<T_scale>::value,
      Eigen::Array<partials_return_type_t<T_scale>, -1, 1>,
      partials_return_type_t<T_scale>>
      T_scale_val;

  using Eigen::Array;
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::VectorXd;

  const size_t N = x.rows();
  const size_t M = x.cols();

  check_positive_finite(function, "Scale vector", sigma);
  check_consistent_size(function, "Vector of dependent variables", y, N);
  check_consistent_size(function, "Weight vector", beta, M);
  if (is_vector<T_scale>::value)
    check_consistent_sizes(function, "Vector of scale parameters", sigma,
                           "Vector of dependent variables", y);
  if (is_vector<T_alpha>::value)
    check_consistent_sizes(function, "Vector of intercepts", alpha,
                           "Vector of dependent variables", y);
  if (size_zero(y, x, beta, sigma))
    return 0;

  if (!include_summand<propto, T_y, T_x, T_alpha, T_beta, T_scale>::value)
    return 0;

  const auto &x_val = value_of_rec(x);
  const auto &beta_val = value_of_rec(beta);
  const auto &alpha_val = value_of_rec(alpha);
  const auto &sigma_val = value_of_rec(sigma);
  const auto &y_val = value_of_rec(y);

  const auto &beta_val_vec = as_column_vector_or_scalar(beta_val);
  const auto &alpha_val_vec = as_column_vector_or_scalar(alpha_val);
  const auto &sigma_val_vec = as_column_vector_or_scalar(sigma_val);
  const auto &y_val_vec = as_column_vector_or_scalar(y_val);

  T_scale_val inv_sigma = 1 / as_array_or_scalar(sigma_val_vec);
  Matrix<T_partials_return, Dynamic, 1> y_minus_mu_over_sigma_mat(N);
  auto y_minus_mu_over_sigma = y_minus_mu_over_sigma_mat.array();

  operands_and_partials<T_y, T_x, T_alpha, T_beta, T_scale> ops_partials(
      y, x, alpha, beta, sigma);

  double y_minus_mu_over_sigma_squared_sum;  // the most efficient way to
  // calculate this depends on
  // template parameters

#ifdef STAN_OPENCL
  const int local_size
      = opencl_kernels::normal_id_glm.make_functor.get_opts().at(
          "LOCAL_SIZE_");
  const int wgs = (N + local_size - 1) / local_size;

  const matrix_cl<double> y_cl = matrix_cl<double>::constant(y_val_vec);
  const matrix_cl<double> x_cl = matrix_cl<double>::constant(x_val);
  matrix_cl<double> beta_cl(beta_val_vec);
  matrix_cl<double> alpha_cl(alpha_val_vec);
  matrix_cl<double> sigma_cl(sigma_val_vec);

  const bool need_mu_derivative
      = !is_constant_all<T_y, T_x, T_beta, T_alpha>::value;
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
  y_minus_mu_over_sigma_squared_sum
      = sum(from_matrix_cl(y_minus_mu_over_sigma_squared_sum_cl));

  if (!is_constant_all<T_y, T_x>::value
      || (!is_constant_all<T_alpha>::value && is_vector<T_alpha>::value)) {
    Matrix<T_partials_return, Dynamic, 1> mu_derivative
        = from_matrix_cl(mu_derivative_cl);
    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_ = -mu_derivative;
    }
    if (!is_constant_all<T_alpha>::value && is_vector<T_alpha>::value) {
      ops_partials.edge3_.partials_ = std::move(mu_derivative);
    }
  }
  if (!is_constant_all<T_alpha>::value && !is_vector<T_alpha>::value) {
    ops_partials.edge3_.partials_[0]
        = sum(from_matrix_cl(mu_derivative_sum_cl));
  }
  if (!is_constant_all<T_x>::value) {
    const matrix_cl<double> beta_transpose_cl(
        beta_cl.buffer(), 1, beta_cl.rows());  // transposition of a vector
                                               // can be done without copying
    ops_partials.edge2_.partials_
        = from_matrix_cl(mu_derivative_cl * beta_transpose_cl);
  }
  if (!is_constant_all<T_beta>::value) {
    const matrix_cl<double> mu_derivative_transpose_cl(
        mu_derivative_cl.buffer(), 1,
        mu_derivative_cl.rows());  // transposition of a vector can be done
                                   // without copying
    ops_partials.edge4_.partials_
        = from_matrix_cl<1, Dynamic>(mu_derivative_transpose_cl * x_cl);
  }
  if (!is_constant_all<T_scale>::value) {
    if (is_vector<T_scale>::value) {
      ops_partials.edge5_.partials_
          = from_matrix_cl<Dynamic, 1>(sigma_derivative_cl);
    } else {
      ops_partials.edge5_.partials_[0]
          = (y_minus_mu_over_sigma_squared_sum - N) * as_scalar(inv_sigma);
    }
  }
#else
  y_minus_mu_over_sigma = x_val * beta_val_vec;
  y_minus_mu_over_sigma
      = (as_array_or_scalar(y_val_vec) - y_minus_mu_over_sigma
         - as_array_or_scalar(alpha_val_vec))
        * inv_sigma;

  if (!(is_constant_all<T_y, T_x, T_beta, T_alpha>::value)) {
    Matrix<T_partials_return, Dynamic, 1> mu_derivative
        = inv_sigma * y_minus_mu_over_sigma;
    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_ = -mu_derivative;
    }
    if (!is_constant_all<T_x>::value) {
      ops_partials.edge2_.partials_
          = (beta_val_vec * mu_derivative.transpose()).transpose();
    }
    if (!is_constant_all<T_beta>::value) {
      ops_partials.edge4_.partials_ = mu_derivative.transpose() * x_val;
    }
    if (!is_constant_all<T_alpha>::value) {
      if (is_vector<T_alpha>::value) {
        ops_partials.edge3_.partials_ = mu_derivative;
      } else {
        ops_partials.edge3_.partials_[0] = sum(mu_derivative);
      }
    }
    if (!is_constant_all<T_scale>::value) {
      if (is_vector<T_scale>::value) {
        Array<T_partials_return, Dynamic, 1> y_minus_mu_over_sigma_squared
            = y_minus_mu_over_sigma * y_minus_mu_over_sigma;
        y_minus_mu_over_sigma_squared_sum
            = sum(y_minus_mu_over_sigma_squared);
        ops_partials.edge5_.partials_
            = (y_minus_mu_over_sigma_squared - 1) * inv_sigma;
      } else {
        y_minus_mu_over_sigma_squared_sum
            = sum(y_minus_mu_over_sigma * y_minus_mu_over_sigma);
        ops_partials.edge5_.partials_[0]
            = (y_minus_mu_over_sigma_squared_sum - N) * as_scalar(inv_sigma);
      }
    }
  } else {
    y_minus_mu_over_sigma_squared_sum
        = sum(y_minus_mu_over_sigma * y_minus_mu_over_sigma);
  }
#endif

  if (!std::isfinite(y_minus_mu_over_sigma_squared_sum)) {
    check_finite(function, "Vector of dependent variables", y);
    check_finite(function, "Weight vector", beta);
    check_finite(function, "Intercept", alpha);
    check_finite(function, "Matrix of independent variables",
                 y_minus_mu_over_sigma_squared_sum);  // if all other checks
                                                      // passed this will only
                                                      // fail if x is not finite
  }

  // Compute log probability.
  T_partials_return logp(0.0);
  if (include_summand<propto>::value)
    logp += NEG_LOG_SQRT_TWO_PI * N;
  if (include_summand<propto, T_scale>::value) {
    if (is_vector<T_scale>::value)
      logp -= sum(log(sigma_val_vec));
    else
      logp -= N * log(as_scalar(sigma_val));
  }
  if (include_summand<propto, T_y, T_x, T_alpha, T_beta, T_scale>::value)
    logp -= 0.5 * y_minus_mu_over_sigma_squared_sum;
  return ops_partials.build(logp);
}

template <typename T_y, typename T_x, typename T_alpha, typename T_beta,
          typename T_scale>
inline return_type_t<T_y, T_x, T_alpha, T_beta, T_scale> normal_id_glm_lpdf(
    const T_y &y, const T_x &x, const T_alpha &alpha, const T_beta &beta,
    const T_scale &sigma) {
  return normal_id_glm_lpdf<false>(y, x, alpha, beta, sigma);
}
}  // namespace math
}  // namespace stan
#endif
