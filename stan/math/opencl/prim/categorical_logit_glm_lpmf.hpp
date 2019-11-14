#ifndef STAN_MATH_OPENCL_PRIM_CATEGORICAL_LOGIT_GLM_LPMF_HPP
#define STAN_MATH_OPENCL_PRIM_CATEGORICAL_LOGIT_GLM_LPMF_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/err/check_bounded.hpp>
#include <stan/math/prim/scal/err/check_consistent_size.hpp>
#include <stan/math/prim/scal/fun/size_zero.hpp>
#include <stan/math/prim/mat/fun/sum.hpp>
#include <Eigen/Core>

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernels/categorical_logit_glm_lpmf.hpp>

namespace stan {
namespace math {

/**
 * Returns the log PMF of the Generalized Linear Model (GLM)
 * with categorical distribution and logit (softmax) link function.
 * This is an overload of the GLM in
 * prim/mar/prob/categorical_logit_glm_lpmf.hpp that is implemented in OpenCL.
 *
 * @tparam T_alpha_scalar type of scalar in the intercept vector
 * @tparam T_beta_scalar type of a scalar in the matrix of weights
 * @param y_cl a scalar or vector of classes. If it is a scalar it will be
 * broadcast - used for all instances. Values should be between 1 and number of
 * classes, including endpoints.
 * @param x_cl design matrix or row vector. This overload does not support
 * broadcasting of a row vector x!
 * @param alpha intercept vector (in log odds)
 * @param beta weight matrix
 * @return log probability or log sum of probabilities
 * @throw std::domain_error x, beta or alpha is infinite or y is not within
 * bounds
 * @throw std::invalid_argument if container sizes mismatch.
 */
template <bool propto, typename T_alpha_scalar, typename T_beta_scalar>
return_type_t<T_alpha_scalar, T_beta_scalar> categorical_logit_glm_lpmf(
    const matrix_cl<int>& y_cl, const matrix_cl<double>& x_cl,
    const Eigen::Matrix<T_alpha_scalar, Eigen::Dynamic, 1>& alpha,
    const Eigen::Matrix<T_beta_scalar, Eigen::Dynamic, Eigen::Dynamic>& beta) {
  using T_partials_return = partials_return_t<T_alpha_scalar, T_beta_scalar>;
  static const char* function = "categorical_logit_glm_lpmf";

  using Eigen::Array;
  using Eigen::Dynamic;
  using Eigen::Matrix;

  const size_t N_instances = x_cl.rows();
  const size_t N_attributes = x_cl.cols();
  const size_t N_classes = beta.cols();

  if (y_cl.size() != 1) {
    check_size_match(function, "x.rows()", N_instances, "y.size()",
                     y_cl.size());
  }
  check_consistent_size(function, "Intercept vector", alpha, N_classes);
  check_size_match(function, "x.cols()", N_attributes, "beta.rows()",
                   beta.rows());

  if (N_instances == 0 || N_classes <= 1) {
    return 0;
  }

  if (!include_summand<propto, T_alpha_scalar, T_beta_scalar>::value) {
    return 0;
  }

  const auto& beta_val = value_of_rec(beta);
  const auto& alpha_val = value_of_rec(alpha);

  const auto& alpha_val_vec = as_column_vector_or_scalar(alpha_val).transpose();

  const int local_size
      = opencl_kernels::categorical_logit_glm.get_option("LOCAL_SIZE_");
  const int wgs = (N_instances + local_size - 1) / local_size;

  const matrix_cl<double> beta_cl(beta_val);
  const matrix_cl<double> alpha_cl(alpha_val);

  matrix_cl<double> x_beta_cl = x_cl * beta_cl;

  bool need_alpha_derivative = !is_constant_all<T_alpha_scalar>::value;
  bool need_beta_derivative = !is_constant_all<T_beta_scalar>::value;

  matrix_cl<double> logp_cl(wgs, 1);
  matrix_cl<double> exp_lin_cl(N_instances, N_classes);
  matrix_cl<double> inv_sum_exp_lin_cl(N_instances, 1);
  matrix_cl<double> neg_softmax_lin_cl(
      need_alpha_derivative || need_beta_derivative ? N_instances : 0,
      N_classes);
  matrix_cl<double> alpha_derivative_cl(need_alpha_derivative ? wgs : 0,
                                        N_classes);

  try {
    opencl_kernels::categorical_logit_glm(
        cl::NDRange(local_size * wgs), cl::NDRange(local_size), logp_cl,
        exp_lin_cl, inv_sum_exp_lin_cl, neg_softmax_lin_cl, alpha_derivative_cl,
        y_cl, x_beta_cl, alpha_cl, N_instances, N_attributes, N_classes,
        y_cl.size() != 1, need_alpha_derivative, need_beta_derivative);
  } catch (const cl::Error& e) {
    check_opencl_error(function, e);
  }
  T_partials_return logp = sum(from_matrix_cl(logp_cl));

  if (!std::isfinite(logp)) {
    check_bounded(function, "categorical outcome out of support",
                  from_matrix_cl(y_cl), 1, N_classes);
    check_finite(function, "Weight vector", beta);
    check_finite(function, "Intercept", alpha);
    check_finite(function, "Matrix of independent variables",
                 from_matrix_cl(x_cl));
  }

  operands_and_partials<Matrix<T_alpha_scalar, Dynamic, 1>,
                        Matrix<T_beta_scalar, Dynamic, Dynamic>>
      ops_partials(alpha, beta);
  if (!is_constant_all<T_alpha_scalar>::value) {
    ops_partials.edge1_.partials_
        = from_matrix_cl(alpha_derivative_cl).colwise().sum();
  }
  if (!is_constant_all<T_beta_scalar>::value && N_attributes != 0) {
    matrix_cl<double> beta_derivative_cl = transpose(x_cl) * neg_softmax_lin_cl;
    matrix_cl<double> temp(N_classes, local_size * N_attributes);
    try {
      opencl_kernels::categorical_logit_glm_beta_derivative(
          cl::NDRange(local_size * N_attributes), cl::NDRange(local_size),
          beta_derivative_cl, temp, y_cl, x_cl, N_instances, N_attributes,
          N_classes, y_cl.size() != 1);
    } catch (const cl::Error& e) {
      check_opencl_error(function, e);
    }
    ops_partials.edge2_.partials_ = from_matrix_cl(beta_derivative_cl);
  }
  return ops_partials.build(logp);
}

template <typename T_alpha_scalar, typename T_beta_scalar>
return_type_t<T_alpha_scalar, T_beta_scalar> categorical_logit_glm_lpmf(
    const matrix_cl<int>& y, const matrix_cl<double>& x,
    const Eigen::Matrix<T_alpha_scalar, Eigen::Dynamic, 1>& alpha,
    const Eigen::Matrix<T_beta_scalar, Eigen::Dynamic, Eigen::Dynamic>& beta) {
  return categorical_logit_glm_lpmf<false>(y, x, alpha, beta);
}

}  // namespace math
}  // namespace stan

#endif
#endif
