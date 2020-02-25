#ifndef STAN_MATH_PRIM_FUN_CHOLESKY_CORR_CONSTRAIN_HPP
#define STAN_MATH_PRIM_FUN_CHOLESKY_CORR_CONSTRAIN_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/log1m.hpp>
#include <stan/math/prim/fun/sqrt.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/fun/corr_constrain.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T, require_t<is_eigen_vector<T>>* = nullptr>
Eigen::Matrix<value_type_t<T>, Eigen::Dynamic, Eigen::Dynamic>
cholesky_corr_constrain(const T& y, int K) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using std::sqrt;
  using T_scalar = value_type_t<T>;
  int k_choose_2 = (K * (K - 1)) / 2;
  check_size_match("cholesky_corr_constrain", "y.size()", y.size(),
                   "k_choose_2", k_choose_2);
  Matrix<T_scalar, Dynamic, 1> z(k_choose_2);
  const Eigen::Ref<const plain_type_t<T>>& y_ref = y;
  for (int i = 0; i < k_choose_2; ++i) {
    z[i] = corr_constrain(y_ref[i]);
  }
  Matrix<T_scalar, Dynamic, Dynamic> x(K, K);
  if (K == 0) {
    return x;
  }
  x.setZero();
  x(0, 0) = 1;
  int k = 0;
  for (int i = 1; i < K; ++i) {
    x(i, 0) = z(k++);
    T_scalar sum_sqs(square(x(i, 0)));
    for (int j = 1; j < i; ++j) {
      x(i, j) = z(k++) * sqrt(1.0 - sum_sqs);
      sum_sqs += square(x(i, j));
    }
    x(i, i) = sqrt(1.0 - sum_sqs);
  }
  return x;
}

// FIXME to match above after debugged
template <typename T, require_t<is_eigen_vector<T>>* = nullptr>
Eigen::Matrix<value_type_t<T>, Eigen::Dynamic, Eigen::Dynamic>
cholesky_corr_constrain(const T& y, int K, value_type_t<T>& lp) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using std::sqrt;
  using T_scalar = value_type_t<T>;
  int k_choose_2 = (K * (K - 1)) / 2;
  check_size_match("cholesky_corr_constrain", "y.size()", y.size(),
                   "k_choose_2", k_choose_2);
  Matrix<T_scalar, Dynamic, 1> z(k_choose_2);
  const Eigen::Ref<const plain_type_t<T>>& y_ref = y;
  for (int i = 0; i < k_choose_2; ++i) {
    z[i] = corr_constrain(y_ref[i], lp);
  }
  Matrix<T_scalar, Dynamic, Dynamic> x(K, K);
  if (K == 0) {
    return x;
  }
  x.setZero();
  x(0, 0) = 1;
  int k = 0;
  for (int i = 1; i < K; ++i) {
    x(i, 0) = z(k++);
    T_scalar sum_sqs = square(x(i, 0));
    for (int j = 1; j < i; ++j) {
      lp += 0.5 * log1m(sum_sqs);
      x(i, j) = z(k++) * sqrt(1.0 - sum_sqs);
      sum_sqs += square(x(i, j));
    }
    x(i, i) = sqrt(1.0 - sum_sqs);
  }
  return x;
}

}  // namespace math
}  // namespace stan
#endif
