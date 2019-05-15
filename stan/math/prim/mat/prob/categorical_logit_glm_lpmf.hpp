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

template<bool propto, typename T_x, typename T_alpha, typename T_beta>
typename return_type<T_x, T_alpha, T_beta>::type
categorical_logit_glm_lpmf(const Eigen::Matrix<int, Eigen::Dynamic, 1>& y,
                           const Eigen::Matrix<T_x, Eigen::Dynamic, Eigen::Dynamic>& x,
                           const T_alpha& alpha, //RowVector<var> or std::vector    .... Matrix?
                           const Eigen::Matrix<T_beta, Eigen::Dynamic, Eigen::Dynamic>& beta) {

  typedef typename stan::partials_return_type<T_x, T_alpha, T_beta>::type T_partials_return;
  static const char* function = "categorical_logit_glm_lpmf";

  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::Array;
  using std::exp;
  using std::log;

  const size_t N_instances = x.rows();
  const size_t N_attributes = x.cols();
  const size_t N_classes = beta.cols();

  check_consistent_size(function, "Vector of dependent variables", y, N_instances);
  check_consistent_size(function, "Intercept vector", alpha, N_classes);
  check_size_match(function, "x.cols()", N_attributes, "beta.rows()", beta.rows());
  check_bounded(function, "categorical outcome out of support", y, 1, N_classes);

  if (size_zero(y, x, beta))
    return 0;

  if (!include_summand<propto, T_x, T_alpha, T_beta>::value)
    return 0;

  const auto& x_val = value_of_rec(x);
  const auto& beta_val = value_of_rec(beta);
  const auto& alpha_val = value_of_rec(alpha);

  const auto& alpha_val_vec = as_column_vector_or_scalar(alpha_val).transpose();

  //TODO logsumexp/softmax overflow
  Array<T_partials_return, Dynamic, Dynamic> lin = (x_val * beta_val).rowwise() + alpha_val_vec; //instances*classes TODO opt multiplication?
  Array<T_partials_return, Dynamic, Dynamic> exp_lin = exp(lin); //instances*classes
  Array<T_partials_return, Dynamic, 1> sum_exp_lin = exp_lin.rowwise().sum(); //instances
  Array<T_partials_return, Dynamic, 1> log_sum_exp_lin = log(sum_exp_lin); //instances

  T_partials_return logp = -log_sum_exp_lin.sum();
  for (int i = 0; i < N_instances; i++) {
    logp += lin(i, y[i] - 1);
  }
  //TODO replace previous block with the following line when we have newer Eigen
  //T_partials_return logp = (lin(Eigen::all,y-1) - log_sum_exp_lin).sum();

  if (!std::isfinite(logp)) {
    check_finite(function, "Weight vector", beta);
    check_finite(function, "Intercept", alpha);
    check_finite(function, "Matrix of independent variables", x);
  }

  Array<T_partials_return, Dynamic, Dynamic> softmax_lin = exp_lin.colwise() / sum_exp_lin; //instances*classes TODO multi instead of div
  Matrix<T_partials_return, 1, Dynamic> sum_softmax_lin = softmax_lin.colwise().sum(); //classes

  // Compute the derivatives.
  operands_and_partials<Matrix<T_x, Dynamic, Dynamic>, T_alpha, Matrix<T_beta, Dynamic, Dynamic>> ops_partials(x, alpha, beta);

  if (!is_constant_struct<T_x>::value) {

    Array<double, Dynamic, Dynamic> beta_y(N_instances, N_attributes);
    for (int i = 0; i < N_instances; i++) {
      beta_y.row(i) = beta_val.col(y[i] - 1);
    }
    ops_partials.edge1_.partials_ = beta_y - (exp_lin.matrix() * beta_val.transpose()).array().colwise() / sum_exp_lin;
    //TODO replace previous block with the following line when we have newer Eigen
    //ops_partials.edge1_.partials_ = beta_val(y-1,all) - (exp_lin*beta.transpose()).colwise() / sum_exp_lin;
  }
  if (!is_constant_struct<T_alpha>::value) {
    ops_partials.edge2_.partials_ = -sum_softmax_lin;
    for (int i = 0; i < N_instances; i++) {
      ops_partials.edge2_.partials_[y[i] - 1] += 1;
    }
  }
  if (!is_constant_struct<T_beta>::value) {

    Matrix<T_partials_return, Dynamic, Dynamic> beta_derivative = -x_val.transpose() * softmax_lin.matrix();
    for (int i = 0; i < N_instances; i++) {
      beta_derivative.col(y[i] - 1) += x_val.row(i);
    }
    ops_partials.edge3_.partials_ = std::move(beta_derivative);
    //TODO replace previous block with the following line when we have newer Eigen
    //ops_partials.edge3_.partials_(Eigen::all, y-1) -= x_val.colwise.sum().transpose();
  }

  return ops_partials.build(logp);

}

template<typename T_x, typename T_alpha, typename T_beta>
typename return_type<T_x, T_alpha, T_beta>::type
categorical_logit_glm_lpmf(const Eigen::Matrix<int, Eigen::Dynamic, 1>& y,
                           const Eigen::Matrix<T_x, Eigen::Dynamic, Eigen::Dynamic>& x,
                           const T_alpha& alpha,
                           const Eigen::Matrix<T_beta, Eigen::Dynamic, Eigen::Dynamic>& beta) {
  return categorical_logit_glm_lpmf<false>(y, x, alpha, beta);
}

}
}

#endif
