#ifndef STAN_MATH_PRIM_MAT_PROB_CATEGORICAL_LOGIT_GLM_LPMF_HPP
#define STAN_MATH_PRIM_MAT_PROB_CATEGORICAL_LOGIT_GLM_LPMF_HPP

#include <stan/math/prim/scal/meta/partials_return_type.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/prim/scal/meta/operands_and_partials.hpp>
#include <stan/math/prim/scal/err/check_consistent_size.hpp>
#include <stan/math/prim/scal/meta/as_array_or_scalar.hpp>
#include <stan/math/prim/mat/meta/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/scal/meta/include_summand.hpp>
#include <stan/math/prim/scal/fun/size_zero.hpp>
#include <Eigen/Core>

namespace stan {
namespace math {

/**
 * Returns the log PMF of the Generalized Linear Model (GLM)
 * with categorical distribution and logit (softmax) link function.
 *
 * @tparam T_x_scalar type of a scalar in the matrix of independent variables
 * (features);
 * @tparam T_alpha type of the intercept(s); This can be Eigen::RowVector or
 * std::vector
 * @tparam T_beta_scalar type of a scalar in the matrix of weights
 * @param y vector of classes; Values should be between 1 and number of classes,
 * including endpoints
 * @param x design matrix
 * @param alpha intercept (in log odds)
 * @param beta weight matrix
 * @return log probability or log sum of probabilities
 * @throw std::domain_error x, beta or alpha is infinite or y is not within
 * bounds
 * @throw std::invalid_argument if container sizes mismatch.
 */
template <bool propto, typename T_x_scalar, typename T_alpha,
          typename T_beta_scalar>
typename return_type<T_x_scalar, T_alpha, T_beta_scalar>::type
categorical_logit_glm_lpmf(
    const Eigen::Matrix<int, Eigen::Dynamic, 1>& y,
    const Eigen::Matrix<T_x_scalar, Eigen::Dynamic, Eigen::Dynamic>& x,
    const T_alpha& alpha,
    const Eigen::Matrix<T_beta_scalar, Eigen::Dynamic, Eigen::Dynamic>& beta) {
  typedef typename stan::partials_return_type<
      T_x_scalar, T_alpha, T_beta_scalar>::type T_partials_return;
  static const char* function = "categorical_logit_glm_lpmf";

  using Eigen::Array;
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using std::exp;
  using std::log;

  const size_t N_instances = x.rows();
  const size_t N_attributes = x.cols();
  const size_t N_classes = beta.cols();

  check_consistent_size(function, "Vector of dependent variables", y,
                        N_instances);
  check_consistent_size(function, "Intercept vector", alpha, N_classes);
  check_size_match(function, "x.cols()", N_attributes, "beta.rows()",
                   beta.rows());
  check_bounded(function, "categorical outcome out of support", y, 1,
                N_classes);

  if (size_zero(y, x, beta))
    return 0;

  if (!include_summand<propto, T_x_scalar, T_alpha, T_beta_scalar>::value)
    return 0;

  const auto& x_val = value_of_rec(x);
  const auto& beta_val = value_of_rec(beta);
  const auto& alpha_val = value_of_rec(alpha);

  const auto& alpha_val_vec = as_column_vector_or_scalar(alpha_val).transpose();

  Array<T_partials_return, Dynamic, Dynamic> lin
      = (x_val * beta_val).rowwise() + alpha_val_vec;
  Array<T_partials_return, Dynamic, 1> lin_max
      = lin.rowwise()
            .maxCoeff();  // This is used to prevent overflow when
                          // calculating softmax/log_sum_exp and similar
                          // expressions
  Array<T_partials_return, Dynamic, Dynamic> exp_lin
      = exp(lin.colwise() - lin_max);
  Array<T_partials_return, Dynamic, 1> inv_sum_exp_lin
      = 1 / exp_lin.rowwise().sum();

  T_partials_return logp = log(inv_sum_exp_lin).sum() - lin_max.sum();
  for (int i = 0; i < N_instances; i++) {
    logp += lin(i, y[i] - 1);
  }
  // TODO(Tadej) maybe we can replace previous block with the following line
  // when we have newer Eigen  T_partials_return logp =
  // lin(Eigen::all,y-1).sum() + log(inv_sum_exp_lin).sum() - lin_max.sum();

  if (!std::isfinite(logp)) {
    check_finite(function, "Weight vector", beta);
    check_finite(function, "Intercept", alpha);
    check_finite(function, "Matrix of independent variables", x);
  }

  // Compute the derivatives.
  operands_and_partials<Matrix<T_x_scalar, Dynamic, Dynamic>, T_alpha,
                        Matrix<T_beta_scalar, Dynamic, Dynamic>>
      ops_partials(x, alpha, beta);

  if (!is_constant_struct<T_x_scalar>::value) {
    Array<double, Dynamic, Dynamic> beta_y(N_instances, N_attributes);
    for (int i = 0; i < N_instances; i++) {
      beta_y.row(i) = beta_val.col(y[i] - 1);
    }
    ops_partials.edge1_.partials_
        = beta_y
          - (exp_lin.matrix() * beta_val.transpose()).array().colwise()
                * inv_sum_exp_lin;
    // TODO(Tadej) maybe we can replace previous block with the following line
    // when we have newer Eigen  ops_partials.edge1_.partials_ = beta_val(y - 1,
    // all) - (exp_lin.matrix() * beta.transpose()).colwise() * inv_sum_exp_lin;
  }
  if (!is_constant_struct<T_alpha>::value
      || !is_constant_struct<T_beta_scalar>::value) {
    Array<T_partials_return, Dynamic, Dynamic> neg_softmax_lin
        = exp_lin.colwise() * -inv_sum_exp_lin;
    if (!is_constant_struct<T_alpha>::value) {
      ops_partials.edge2_.partials_ = neg_softmax_lin.colwise().sum();
      for (int i = 0; i < N_instances; i++) {
        ops_partials.edge2_.partials_[y[i] - 1] += 1;
      }
    }
    if (!is_constant_struct<T_beta_scalar>::value) {
      Matrix<T_partials_return, Dynamic, Dynamic> beta_derivative
          = x_val.transpose() * neg_softmax_lin.matrix();

      for (int i = 0; i < N_instances; i++) {
        beta_derivative.col(y[i] - 1) += x_val.row(i);
      }
      // TODO(Tadej) maybe we can replace previous loop with the following line
      // when we have newer Eigen  ops_partials.edge3_.partials_(Eigen::all, y -
      // 1) += x_val.colwise.sum().transpose();

      ops_partials.edge3_.partials_ = std::move(beta_derivative);
    }
  }
  return ops_partials.build(logp);
}

template <typename T_x, typename T_alpha, typename T_beta>
typename return_type<T_x, T_alpha, T_beta>::type categorical_logit_glm_lpmf(
    const Eigen::Matrix<int, Eigen::Dynamic, 1>& y,
    const Eigen::Matrix<T_x, Eigen::Dynamic, Eigen::Dynamic>& x,
    const T_alpha& alpha,
    const Eigen::Matrix<T_beta, Eigen::Dynamic, Eigen::Dynamic>& beta) {
  return categorical_logit_glm_lpmf<false>(y, x, alpha, beta);
}

}  // namespace math
}  // namespace stan

#endif
