#ifndef STAN_MATH_OPENCL_PRIM_CATEGORICAL_LOGIT_GLM_LPMF_HPP
#define STAN_MATH_OPENCL_PRIM_CATEGORICAL_LOGIT_GLM_LPMF_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/prim/size.hpp>
#include <stan/math/opencl/rev/arena_matrix_cl.hpp>
#include <stan/math/opencl/rev/operands_and_partials.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/copy.hpp>
#include <stan/math/opencl/prim/multiply.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/kernels/categorical_logit_glm_lpmf.hpp>

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/eval.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/** \ingroup opencl
 * Returns the log PMF of the Generalized Linear Model (GLM)
 * with categorical distribution and logit (softmax) link function.
 * This is an overload of the GLM in
 * prim/prob/categorical_logit_glm_lpmf.hpp that is implemented in OpenCL.
 *
 * @tparam T_alpha type of the intercept vector
 * @tparam T_beta type of the matrix of weights
 * @param y a scalar or vector of classes. If it is a scalar it will be
 * broadcast - used for all instances. Values should be between 1 and number of
 * classes, including endpoints.
 * @param x design matrix on OpenCL device. This overload does not support
 * broadcasting of a row vector x!
 * @param alpha intercept vector (in log odds)
 * @param beta weight matrix
 * @return log probability or log sum of probabilities
 * @throw std::domain_error x, beta or alpha is infinite or y is not within
 * bounds
 * @throw std::invalid_argument if container sizes mismatch.
 */
template <bool propto, typename T_y, typename T_x, typename T_alpha,
          typename T_beta,
          require_all_prim_or_rev_kernel_expression_t<T_y, T_x, T_alpha,
                                                      T_beta>* = nullptr>
return_type_t<T_x, T_alpha, T_beta> categorical_logit_glm_lpmf(
    const T_y& y, const T_x& x, const T_alpha& alpha, const T_beta& beta) {
  using T_partials_return = partials_return_t<T_x, T_alpha, T_beta>;
  constexpr bool is_y_vector = !is_stan_scalar<T_y>::value;
  using Eigen::Array;
  using Eigen::Dynamic;
  using Eigen::Matrix;

  const size_t N_instances = x.rows();
  const size_t N_attributes = x.cols();
  const size_t N_classes = beta.cols();

  static const char* function = "categorical_logit_glm_lpmf";
  if (is_y_vector) {
    check_size_match(function, "Rows of ", "x", N_instances, "size of ", "y",
                     math::size(y));
  }
  check_size_match(function, "Columns of ", "beta", N_classes, "size of ",
                   "alpha", math::size(alpha));
  check_size_match(function, "Columns of ", "x", N_attributes, "Rows of",
                   "beta", beta.rows());

  if (N_instances == 0 || N_classes <= 1) {
    return 0;
  }
  if (!include_summand<propto, T_x, T_alpha, T_beta>::value) {
    return 0;
  }

  const auto& y_val = eval(value_of(y));
  const auto& x_val = eval(value_of(x));
  const auto& alpha_val = eval(value_of(alpha));
  const auto& beta_val = eval(value_of(beta));

  const auto& y_val_cl = to_matrix_cl(y_val);

  matrix_cl<double> x_beta_cl = x_val * beta_val;

  const int local_size
      = opencl_kernels::categorical_logit_glm.get_option("LOCAL_SIZE_");
  const int wgs = (N_instances + local_size - 1) / local_size;

  bool need_alpha_derivative = !is_constant_all<T_alpha>::value;
  bool need_beta_derivative = !is_constant_all<T_beta>::value;

  matrix_cl<double> logp_cl(wgs, 1);
  matrix_cl<double> exp_lin_cl(N_instances, N_classes);
  matrix_cl<double> inv_sum_exp_lin_cl(N_instances, 1);
  matrix_cl<double> neg_softmax_lin_cl(
      need_alpha_derivative || need_beta_derivative ? N_instances : 0,
      N_classes);
  matrix_cl<double> alpha_derivative_cl(N_classes,
                                        need_alpha_derivative ? wgs : 0);

  try {
    opencl_kernels::categorical_logit_glm(
        cl::NDRange(local_size * wgs), cl::NDRange(local_size), logp_cl,
        exp_lin_cl, inv_sum_exp_lin_cl, neg_softmax_lin_cl, alpha_derivative_cl,
        y_val_cl, x_beta_cl, alpha_val, N_instances, N_attributes, N_classes,
        is_y_vector, need_alpha_derivative, need_beta_derivative);
  } catch (const cl::Error& e) {
    check_opencl_error(function, e);
  }
  T_partials_return logp = sum(from_matrix_cl(logp_cl));

  if (!std::isfinite(logp)) {
    results(check_cl(function, "Vector of dependent variables", y_val,
                     "between 0 and cols of beta"),
            check_cl(function, "Intercept", alpha_val, "finite"))
        = expressions(y_val >= 0 && y_val <= static_cast<int>(N_classes),
                      isfinite(alpha_val));
    check_cl(function, "Design matrix", x_val, "finite") = isfinite(x_val);
    check_cl(function, "Weight vector", beta_val, "finite")
        = isfinite(beta_val);
  }

  auto ops_partials = make_partials_propagator(x, alpha, beta);
  if (!is_constant_all<T_x>::value) {
    if (is_y_vector) {
      partials<0>(ops_partials)
          = indexing(beta_val, col_index(x.rows(), x.cols()),
                     rowwise_broadcast(forward_as<matrix_cl<int>>(y_val) - 1))
            - elt_multiply(exp_lin_cl * transpose(beta_val),
                           rowwise_broadcast(inv_sum_exp_lin_cl));
    } else {
      partials<0>(ops_partials)
          = indexing(beta_val, col_index(x.rows(), x.cols()),
                     forward_as<int>(y_val) - 1)
            - elt_multiply(exp_lin_cl * transpose(beta_val),
                           rowwise_broadcast(inv_sum_exp_lin_cl));
    }
  }
  if (!is_constant_all<T_alpha>::value) {
    if (wgs == 1) {
      partials<1>(ops_partials) = std::move(alpha_derivative_cl);
    } else {
      partials<1>(ops_partials) = rowwise_sum(alpha_derivative_cl);
    }
  }
  if (!is_constant_all<T_beta>::value && N_attributes != 0) {
    partials<2>(ops_partials) = transpose(x_val) * neg_softmax_lin_cl;
    matrix_cl<double> temp(N_classes, local_size * N_attributes);
    try {
      opencl_kernels::categorical_logit_glm_beta_derivative(
          cl::NDRange(local_size * N_attributes), cl::NDRange(local_size),
          forward_as<arena_matrix_cl<double>>(partials<2>(ops_partials)), temp,
          y_val_cl, x_val, N_instances, N_attributes, N_classes, is_y_vector);
    } catch (const cl::Error& e) {
      check_opencl_error(function, e);
    }
  }
  return ops_partials.build(logp);
}

}  // namespace math
}  // namespace stan

#endif
#endif
