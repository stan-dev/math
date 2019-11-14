#ifndef STAN_MATH_PRIM_MAT_PROB_CATEGORICAL_LOGIT_GLM_LPMF_HPP
#define STAN_MATH_PRIM_MAT_PROB_CATEGORICAL_LOGIT_GLM_LPMF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/err/check_bounded.hpp>
#include <stan/math/prim/scal/err/check_consistent_size.hpp>
#include <stan/math/prim/scal/fun/size_zero.hpp>
#include <Eigen/Core>

namespace stan {
namespace math {

/**
 * Returns the log PMF of the Generalized Linear Model (GLM)
 * with categorical distribution and logit (softmax) link function.
 *
 * @tparam T_y type of classes. It can be either `std::vector<int>` or `int`.
 * @tparam T_x_scalar type of a scalar in the matrix of independent variables
 * (features)
 * @tparam T_x_rows compile-time number of rows of `x`. It can be either
 * `Eigen::Dynamic` or 1.
 * @tparam T_alpha_scalar type of scalar in the intercept vector
 * @tparam T_beta_scalar type of a scalar in the matrix of weights
 * @param y a scalar or vector of classes. If it is a scalar it will be
 * broadcast - used for all instances. Values should be between 1 and number of
 * classes, including endpoints.
 * @param x design matrix or row vector. If it is a row vector it will be
 * broadcast - used for all instances.
 * @param alpha intercept vector (in log odds)
 * @param beta weight matrix
 * @return log probability or log sum of probabilities
 * @throw std::domain_error x, beta or alpha is infinite or y is not within
 * bounds
 * @throw std::invalid_argument if container sizes mismatch.
 */
template <bool propto, typename T_y, typename T_x_scalar, int T_x_rows,
          typename T_alpha_scalar, typename T_beta_scalar>
return_type_t<T_x_scalar, T_alpha_scalar, T_beta_scalar>
categorical_logit_glm_lpmf(
    const T_y& y, const Eigen::Matrix<T_x_scalar, T_x_rows, Eigen::Dynamic>& x,
    const Eigen::Matrix<T_alpha_scalar, Eigen::Dynamic, 1>& alpha,
    const Eigen::Matrix<T_beta_scalar, Eigen::Dynamic, Eigen::Dynamic>& beta) {
  using T_partials_return
      = partials_return_t<T_x_scalar, T_alpha_scalar, T_beta_scalar>;
  static const char* function = "categorical_logit_glm_lpmf";

  using Eigen::Array;
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using std::exp;
  using std::log;

  const size_t N_instances = T_x_rows == 1 ? length(y) : x.rows();
  const size_t N_attributes = x.cols();
  const size_t N_classes = beta.cols();

  if (is_vector<T_y>::value && T_x_rows != 1) {
    check_consistent_size(function, "Vector of dependent variables", y,
                          N_instances);
  }
  check_consistent_size(function, "Intercept vector", alpha, N_classes);
  check_size_match(function, "x.cols()", N_attributes, "beta.rows()",
                   beta.rows());
  check_bounded(function, "categorical outcome out of support", y, 1,
                N_classes);

  if (size_zero(y) || N_classes == 1) {
    return 0;
  }

  if (!include_summand<propto, T_x_scalar, T_alpha_scalar,
                       T_beta_scalar>::value) {
    return 0;
  }

  const auto& x_val = value_of_rec(x);
  const auto& beta_val = value_of_rec(beta);
  const auto& alpha_val = value_of_rec(alpha);

  const auto& alpha_val_vec = as_column_vector_or_scalar(alpha_val).transpose();

  Array<T_partials_return, T_x_rows, Dynamic> lin
      = (x_val * beta_val).rowwise() + alpha_val_vec;
  Array<T_partials_return, T_x_rows, 1> lin_max
      = lin.rowwise().maxCoeff();  // This is used to prevent overflow when
                                   // calculating softmax/log_sum_exp and
                                   // similar expressions
  Array<T_partials_return, T_x_rows, Dynamic> exp_lin
      = exp(lin.colwise() - lin_max);
  Array<T_partials_return, T_x_rows, 1> inv_sum_exp_lin
      = 1 / exp_lin.rowwise().sum();

  T_partials_return logp = log(inv_sum_exp_lin).sum() - lin_max.sum();
  if (T_x_rows == 1) {
    logp *= N_instances;
  }
  scalar_seq_view<T_y> y_seq(y);
  for (int i = 0; i < N_instances; i++) {
    if (T_x_rows == 1) {
      logp += lin(0, y_seq[i] - 1);
    } else {
      logp += lin(i, y_seq[i] - 1);
    }
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
  operands_and_partials<Matrix<T_x_scalar, T_x_rows, Dynamic>,
                        Matrix<T_alpha_scalar, Dynamic, 1>,
                        Matrix<T_beta_scalar, Dynamic, Dynamic>>
      ops_partials(x, alpha, beta);

  if (!is_constant_all<T_x_scalar>::value) {
    if (T_x_rows == 1) {
      Array<double, 1, Dynamic> beta_y = beta_val.col(y_seq[0] - 1);
      for (int i = 1; i < N_instances; i++) {
        beta_y += beta_val.col(y_seq[i] - 1).array();
      }
      ops_partials.edge1_.partials_
          = beta_y
            - (exp_lin.matrix() * beta_val.transpose()).array().colwise()
                  * inv_sum_exp_lin * N_instances;
    } else {
      Array<double, Dynamic, Dynamic> beta_y(N_instances, N_attributes);
      for (int i = 0; i < N_instances; i++) {
        beta_y.row(i) = beta_val.col(y_seq[i] - 1);
      }
      ops_partials.edge1_.partials_
          = beta_y
            - (exp_lin.matrix() * beta_val.transpose()).array().colwise()
                  * inv_sum_exp_lin;
      // TODO(Tadej) maybe we can replace previous block with the following line
      // when we have newer Eigen  ops_partials.edge1_.partials_ = beta_val(y -
      // 1, all) - (exp_lin.matrix() * beta.transpose()).colwise() *
      // inv_sum_exp_lin;
    }
  }
  if (!is_constant_all<T_alpha_scalar, T_beta_scalar>::value) {
    Array<T_partials_return, T_x_rows, Dynamic> neg_softmax_lin
        = exp_lin.colwise() * -inv_sum_exp_lin;
    if (!is_constant_all<T_alpha_scalar>::value) {
      if (T_x_rows == 1) {
        ops_partials.edge2_.partials_
            = neg_softmax_lin.colwise().sum() * N_instances;
      } else {
        ops_partials.edge2_.partials_ = neg_softmax_lin.colwise().sum();
      }
      for (int i = 0; i < N_instances; i++) {
        ops_partials.edge2_.partials_[y_seq[i] - 1] += 1;
      }
    }
    if (!is_constant_all<T_beta_scalar>::value) {
      Matrix<T_partials_return, Dynamic, Dynamic> beta_derivative;
      if (T_x_rows == 1) {
        beta_derivative
            = x_val.transpose() * neg_softmax_lin.matrix() * N_instances;
      } else {
        beta_derivative = x_val.transpose() * neg_softmax_lin.matrix();
      }

      for (int i = 0; i < N_instances; i++) {
        if (T_x_rows == 1) {
          beta_derivative.col(y_seq[i] - 1) += x_val;
        } else {
          beta_derivative.col(y_seq[i] - 1) += x_val.row(i);
        }
      }
      // TODO(Tadej) maybe we can replace previous loop with the following line
      // when we have newer Eigen  ops_partials.edge3_.partials_(Eigen::all, y -
      // 1) += x_val.colwise.sum().transpose();

      ops_partials.edge3_.partials_ = std::move(beta_derivative);
    }
  }
  return ops_partials.build(logp);
}

template <typename T_y, typename T_x_scalar, int T_x_rows,
          typename T_alpha_scalar, typename T_beta_scalar>
return_type_t<T_x_scalar, T_alpha_scalar, T_beta_scalar>
categorical_logit_glm_lpmf(
    const T_y& y, const Eigen::Matrix<T_x_scalar, T_x_rows, Eigen::Dynamic>& x,
    const Eigen::Matrix<T_alpha_scalar, Eigen::Dynamic, 1>& alpha,
    const Eigen::Matrix<T_beta_scalar, Eigen::Dynamic, Eigen::Dynamic>& beta) {
  return categorical_logit_glm_lpmf<false>(y, x, alpha, beta);
}

}  // namespace math
}  // namespace stan

#endif
