#ifndef STAN_MATH_OPENCL_PRIM_BERNOULLI_LOGIT_GLM_LPMF_HPP
#define STAN_MATH_OPENCL_PRIM_BERNOULLI_LOGIT_GLM_LPMF_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_bounded.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/mat/fun/value_of_rec.hpp>
#include <stan/math/prim/mat/fun/sum.hpp>
#include <stan/math/prim/arr/fun/value_of_rec.hpp>
#include <stan/math/prim/scal/fun/size_zero.hpp>

#include <stan/math/prim/mat/fun/value_of.hpp>
#include <stan/math/prim/arr/fun/value_of.hpp>
#include <stan/math/opencl/kernels/bernoulli_logit_glm_lpmf.hpp>

#include <cmath>

namespace stan {
namespace math {

/**
 * Returns the log PMF of the Generalized Linear Model (GLM)
 * with Bernoulli distribution and logit link function.
 * This is an overload of the GLM in prim/mat/prob/bernoulli_logit_glm_lpmf.hpp
 * that is implemented in OpenCL.
 * @tparam T_alpha type of the intercept(s);
 * this can be a vector (of the same length as y) of intercepts or a single
 * value (for models with constant intercept);
 * @tparam T_beta type of the weight vector;
 * this can also be a single value;
 * @param y_cl binary vector parameter on OpenCL device
 * @param x_cl design matrix on OpenCL device
 * @param alpha intercept (in log odds)
 * @param beta weight vector
 * @return log probability or log sum of probabilities
 * @throw std::domain_error if x, beta or alpha is infinite.
 * @throw std::domain_error if y is not binary.
 * @throw std::invalid_argument if container sizes mismatch.
 */

template <bool propto, typename T_alpha, typename T_beta>
return_type_t<T_alpha, T_beta> bernoulli_logit_glm_lpmf(
    const matrix_cl<int> &y_cl, const matrix_cl<double> &x_cl,
    const T_alpha &alpha, const T_beta &beta) {
  static const char *function = "bernoulli_logit_glm_lpmf(OpenCL)";
  using T_partials_return = partials_return_t<T_alpha, T_beta>;

  using Eigen::Dynamic;
  using Eigen::Matrix;

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
  const auto &beta_val = value_of_rec(beta);
  const auto &alpha_val = value_of_rec(alpha);

  const auto &beta_val_vec = as_column_vector_or_scalar(beta_val);
  const auto &alpha_val_vec = as_column_vector_or_scalar(alpha_val);

  const int local_size
      = opencl_kernels::bernoulli_logit_glm.get_option("LOCAL_SIZE_");
  const int wgs = (N + local_size - 1) / local_size;

  matrix_cl<double> beta_cl(beta_val_vec);
  matrix_cl<double> alpha_cl(alpha_val_vec);

  matrix_cl<double> logp_cl(wgs, 1);
  const bool need_theta_derivative = !is_constant_all<T_beta, T_alpha>::value;
  matrix_cl<double> theta_derivative_cl(need_theta_derivative ? N : 0, 1);
  const bool need_theta_derivative_sum
      = need_theta_derivative && !is_vector<T_alpha>::value;
  matrix_cl<double> theta_derivative_sum_cl(need_theta_derivative_sum ? wgs : 0,
                                            1);

  try {
    opencl_kernels::bernoulli_logit_glm(
        cl::NDRange(local_size * wgs), cl::NDRange(local_size), logp_cl,
        theta_derivative_cl, theta_derivative_sum_cl, y_cl, x_cl, alpha_cl,
        beta_cl, N, M, length(alpha) != 1, need_theta_derivative,
        need_theta_derivative_sum);
  } catch (const cl::Error &e) {
    check_opencl_error(function, e);
  }

  Eigen::VectorXd logp_partial_sum(wgs);
  logp_partial_sum = from_matrix_cl(logp_cl);
  logp += sum(logp_partial_sum);

  if (!std::isfinite(logp)) {
    check_bounded(function, "Vector of dependent variables",
                  from_matrix_cl(y_cl), 0, 1);
    check_finite(function, "Weight vector", beta);
    check_finite(function, "Intercept", alpha);
    check_finite(function, "Matrix of independent variables",
                 from_matrix_cl(x_cl));
  }

  operands_and_partials<T_alpha, T_beta> ops_partials(alpha, beta);
  // Compute the necessary derivatives.
  if (!is_constant_all<T_alpha>::value) {
    if (is_vector<T_alpha>::value) {
      ops_partials.edge1_.partials_
          = std::move(from_matrix_cl<Dynamic, 1>(theta_derivative_cl));
    } else {
      ops_partials.edge1_.partials_[0]
          = sum(from_matrix_cl(theta_derivative_sum_cl));
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
  return ops_partials.build(logp);
}

template <typename T_alpha, typename T_beta>
inline return_type_t<T_beta, T_alpha> bernoulli_logit_glm_lpmf(
    const matrix_cl<int> &y, const matrix_cl<double> &x, const T_alpha &alpha,
    const T_beta &beta) {
  return bernoulli_logit_glm_lpmf<false>(y, x, alpha, beta);
}

}  // namespace math
}  // namespace stan

#endif
#endif
