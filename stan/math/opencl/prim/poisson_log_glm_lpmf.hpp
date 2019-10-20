#ifndef STAN_MATH_OPENCL_PRIM_POISSON_LOG_GLM_LPMF_HPP
#define STAN_MATH_OPENCL_PRIM_POISSON_LOG_GLM_LPMF_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_nonnegative.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/fun/lgamma.hpp>
#include <stan/math/prim/scal/fun/size_zero.hpp>
#include <stan/math/prim/mat/fun/value_of_rec.hpp>
#include <stan/math/prim/arr/fun/value_of_rec.hpp>
#include <stan/math/prim/mat/fun/value_of.hpp>
#include <stan/math/prim/arr/fun/value_of.hpp>
#include <stan/math/prim/mat/fun/sum.hpp>
#include <stan/math/opencl/kernels/poisson_log_glm_lpmf.hpp>
#include <cmath>
#include <limits>

namespace stan {
namespace math {

/**
 * Returns the log PMF of the Generalized Linear Model (GLM)
 * with Poisson distribution and log link function.
 * This is an overload of the GLM in prim/mat/prob/poisson_log_glm_lpmf.hpp
 * that is implemented in OpenCL.
 * @tparam T_alpha type of the intercept(s);
 * this can be a vector (of the same length as y) of intercepts or a single
 * value (for models with constant intercept);
 * @tparam T_beta type of the weight vector;
 * this can also be a single value;
 * @param y_cl positive integer vector parameter on OpenCL device
 * @param x_cl design matrix on OpenCL device
 * @param alpha intercept (in log odds)
 * @param beta weight vector
 * @return log probability or log sum of probabilities
 * @throw std::domain_error if x, beta or alpha is infinite.
 * @throw std::domain_error if y is negative.
 * @throw std::invalid_argument if container sizes mismatch.
 */
template <bool propto, typename T_alpha, typename T_beta>
return_type_t<T_alpha, T_beta> poisson_log_glm_lpmf(
    const matrix_cl<int>& y_cl, const matrix_cl<double>& x_cl,
    const T_alpha& alpha, const T_beta& beta) {
  static const char* function = "poisson_log_glm_lpmf(OpenCL)";
  using T_partials_return = partials_return_t<T_alpha, T_beta>;

  using Eigen::Dynamic;
  using Eigen::Matrix;
  using std::exp;

  const size_t N = x_cl.rows();
  const size_t M = x_cl.cols();

  check_size_match(function, "Rows of ", "x_cl", N, "rows of ", "y_cl",
                   y_cl.rows());
  check_consistent_size(function, "Weight vector", beta, M);
  if (is_vector<T_alpha>::value) {
    check_size_match(function, "Rows of ", "y_cl", N, "size of ", "alpha",
                     length(alpha));
  }
  if (N == 0) {
    return 0;
  }

  if (!include_summand<propto, T_alpha, T_beta>::value) {
    return 0;
  }

  T_partials_return logp(0);

  const auto& beta_val = value_of_rec(beta);
  const auto& alpha_val = value_of_rec(alpha);

  const auto& beta_val_vec = as_column_vector_or_scalar(beta_val);
  const auto& alpha_val_vec = as_column_vector_or_scalar(alpha_val);

  const int local_size
      = opencl_kernels::poisson_log_glm.get_option("LOCAL_SIZE_");
  const int wgs = (N + local_size - 1) / local_size;

  matrix_cl<double> beta_cl(beta_val_vec);
  matrix_cl<double> alpha_cl(alpha_val_vec);

  matrix_cl<double> theta_derivative_cl(N, 1);
  matrix_cl<double> theta_derivative_sum_cl(wgs, 1);
  const bool need_logp1 = include_summand<propto>::value;
  const bool need_logp2 = include_summand<propto, T_partials_return>::value;
  matrix_cl<double> logp_cl((need_logp1 || need_logp2) ? wgs : 0, 1);

  try {
    opencl_kernels::poisson_log_glm(
        cl::NDRange(local_size * wgs), cl::NDRange(local_size),
        theta_derivative_cl, theta_derivative_sum_cl, logp_cl, y_cl, x_cl,
        alpha_cl, beta_cl, N, M, length(alpha) != 1, need_logp1, need_logp2);
  } catch (const cl::Error& e) {
    check_opencl_error(function, e);
  }
  Matrix<T_partials_return, Dynamic, 1> theta_derivative_partial_sum(wgs);
  theta_derivative_partial_sum = from_matrix_cl(theta_derivative_sum_cl);
  double theta_derivative_sum = sum(theta_derivative_partial_sum);
  if (need_logp1 || need_logp2) {
    Eigen::VectorXd logp_partial_sum(wgs);
    logp_partial_sum = from_matrix_cl(logp_cl);
    logp += sum(logp_partial_sum);
  }
  if (!std::isfinite(theta_derivative_sum)) {
    check_nonnegative(function, "Vector of dependent variables",
                      from_matrix_cl(y_cl));
    check_finite(function, "Weight vector", beta);
    check_finite(function, "Intercept", alpha);
    check_finite(function, "Matrix of independent variables",
                 from_matrix_cl(x_cl));
  }

  operands_and_partials<T_alpha, T_beta> ops_partials(alpha, beta);
  // Compute the necessary derivatives.
  if (!is_constant_all<T_alpha>::value) {
    if (is_vector<T_alpha>::value)
      ops_partials.edge1_.partials_
          = from_matrix_cl<Dynamic, 1>(theta_derivative_cl);
    else
      ops_partials.edge1_.partials_[0] = theta_derivative_sum;
  }
  if (!is_constant_all<T_beta>::value) {
    matrix_cl<double> theta_derivative_transpose_cl(
        theta_derivative_cl.buffer(), 1,
        theta_derivative_cl
            .rows());  // transposition of a vector can be done without copying
    ops_partials.edge2_.partials_
        = from_matrix_cl<1, Dynamic>(theta_derivative_transpose_cl * x_cl);
  }
  return ops_partials.build(logp);
}

template <typename T_alpha, typename T_beta>
inline return_type_t<T_alpha, T_beta> poisson_log_glm_lpmf(
    const matrix_cl<int>& y, const matrix_cl<double>& x, const T_alpha& alpha,
    const T_beta& beta) {
  return poisson_log_glm_lpmf<false>(y, x, alpha, beta);
}
}  // namespace math
}  // namespace stan

#endif
#endif
