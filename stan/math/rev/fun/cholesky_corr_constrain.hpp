#ifndef STAN_MATH_REV_FUN_CHOLESKY_CORR_CONSTRAIN_HPP
#define STAN_MATH_REV_FUN_CHOLESKY_CORR_CONSTRAIN_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/log1m.hpp>
#include <stan/math/prim/fun/sqrt.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/fun/corr_constrain.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the Cholesky factor of the correlation matrix of the sepcified
 * size read from the unconstrained vector `y`. A total of K choose 2
 * elements are required to build a K by K Cholesky factor.
 *
 * @tparam T type of input vector (must be a `var_value<S>` where `S`
 *  inherits from EigenBase)
 * @param y Vector of unconstrained values
 * @param K number of rows
 * @return Cholesky factor of correlation matrix
 */
template <typename T, require_var_vector_t<T>* = nullptr>
var_value<Eigen::MatrixXd> cholesky_corr_constrain(const T& y, int K) {
  using std::sqrt;
  int k_choose_2 = (K * (K - 1)) / 2;
  check_size_match("cholesky_corr_constrain", "y.size()", y.size(),
                   "k_choose_2", k_choose_2);

  if (K == 0) {
    return Eigen::MatrixXd();
  }

  var_value<Eigen::VectorXd> z = corr_constrain(y);
  arena_t<Eigen::MatrixXd> x_val = Eigen::MatrixXd::Zero(K, K);

  x_val.coeffRef(0, 0) = 1;
  int k = 0;
  for (int i = 1; i < K; ++i) {
    x_val.coeffRef(i, 0) = z.val().coeff(k);
    k++;
    double sum_sqs = square(x_val.coeff(i, 0));
    for (int j = 1; j < i; ++j) {
      x_val.coeffRef(i, j) = z.val().coeff(k) * sqrt(1.0 - sum_sqs);
      k++;
      sum_sqs += square(x_val.coeff(i, j));
    }
    x_val.coeffRef(i, i) = sqrt(1.0 - sum_sqs);
  }

  var_value<Eigen::MatrixXd> x = x_val;

  reverse_pass_callback([z, K, x]() mutable {
    size_t k = z.size();
    for (int i = K - 1; i > 0; --i) {
      double sum_sqs_val = 1.0 - square(x.val().coeffRef(i, i));
      double sum_sqs_adj = -0.5 * x.adj().coeff(i, i) / x.val().coeff(i, i);
      for (int j = i - 1; j > 0; --j) {
        x.adj().coeffRef(i, j) += 2 * sum_sqs_adj * x.val().coeff(i, j);
        sum_sqs_val -= square(x.val().coeff(i, j));
        k--;
        sum_sqs_adj += -0.5 * x.adj().coeffRef(i, j) * z.val().coeff(k)
                       / sqrt(1.0 - sum_sqs_val);
        z.adj().coeffRef(k) += x.adj().coeffRef(i, j) * sqrt(1.0 - sum_sqs_val);
      }
      x.adj().coeffRef(i, 0) += 2 * sum_sqs_adj * x.val().coeff(i, 0);
      k--;
      z.adj().coeffRef(k) += x.adj().coeffRef(i, 0);
    }
  });

  return x;
}

/**
 * Return the Cholesky factor of the correlation matrix of the sepcified
 * size read from the unconstrained vector `y`. A total of K choose 2
 * elements are required to build a K by K Cholesky factor.
 *
 * @tparam T type of input vector (must be a `var_value<S>` where `S`
 *  inherits from EigenBase)
 * @param y Vector of unconstrained values
 * @param K number of rows
 * @param[out] lp Log density that is incremented with the log Jacobian
 * @return Cholesky factor of correlation matrix
 */
template <typename T, require_var_vector_t<T>* = nullptr>
var_value<Eigen::MatrixXd> cholesky_corr_constrain(const T& y, int K,
                                                   scalar_type_t<T>& lp) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using std::sqrt;
  int k_choose_2 = (K * (K - 1)) / 2;
  check_size_match("cholesky_corr_constrain", "y.size()", y.size(),
                   "k_choose_2", k_choose_2);

  if (K == 0) {
    return Eigen::MatrixXd();
  }

  var_value<Eigen::VectorXd> z = corr_constrain(y, lp);
  arena_t<Eigen::MatrixXd> x_val = Eigen::MatrixXd::Zero(K, K);

  x_val.coeffRef(0, 0) = 1;
  int k = 0;
  double lp_val = 0.0;
  for (int i = 1; i < K; ++i) {
    x_val.coeffRef(i, 0) = z.val().coeff(k);
    k++;
    double sum_sqs = square(x_val.coeff(i, 0));
    for (int j = 1; j < i; ++j) {
      lp_val += 0.5 * log1m(sum_sqs);
      x_val.coeffRef(i, j) = z.val().coeff(k) * sqrt(1.0 - sum_sqs);
      k++;
      sum_sqs += square(x_val.coeff(i, j));
    }
    x_val.coeffRef(i, i) = sqrt(1.0 - sum_sqs);
  }

  lp += lp_val;
  var_value<Eigen::MatrixXd> x = x_val;

  reverse_pass_callback([z, K, x, lp]() mutable {
    size_t k = z.size();
    for (int i = K - 1; i > 0; --i) {
      double sum_sqs_val = 1.0 - square(x.val().coeffRef(i, i));
      double sum_sqs_adj = -0.5 * x.adj().coeff(i, i) / x.val().coeff(i, i);
      for (int j = i - 1; j > 0; --j) {
        x.adj().coeffRef(i, j) += 2 * sum_sqs_adj * x.val().coeff(i, j);
        sum_sqs_val -= square(x.val().coeff(i, j));

        k--;
        sum_sqs_adj += -0.5 * x.adj().coeffRef(i, j) * z.val().coeff(k)
                       / sqrt(1.0 - sum_sqs_val);
        z.adj().coeffRef(k) += x.adj().coeffRef(i, j) * sqrt(1.0 - sum_sqs_val);
        sum_sqs_adj -= 0.5 * lp.adj() / (1 - sum_sqs_val);
      }
      x.adj().coeffRef(i, 0) += 2 * sum_sqs_adj * x.val().coeff(i, 0);
      k--;
      z.adj().coeffRef(k) += x.adj().coeffRef(i, 0);
    }
  });

  return x;
}

}  // namespace math
}  // namespace stan
#endif
