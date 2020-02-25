#ifndef STAN_MATH_PRIM_FUN_CHOLESKY_CORR_FREE_HPP
#define STAN_MATH_PRIM_FUN_CHOLESKY_CORR_FREE_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/corr_free.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T, require_eigen_t<T>* = nullptr>
Eigen::Matrix<value_type_t<T>, Eigen::Dynamic, 1> cholesky_corr_free(
    const T& x) {
  using Eigen::Dynamic;
  using Eigen::Matrix;

  check_square("cholesky_corr_free", "x", x);
  // should validate lower-triangular, unit lengths

  const Eigen::Ref<const plain_type_t<T>> x_ref = x;
  int K = (x.rows() * (x.rows() - 1)) / 2;
  Matrix<value_type_t<T>, Dynamic, 1> z(K);
  int k = 0;
  for (int i = 1; i < x.rows(); ++i) {
    z(k++) = corr_free(x_ref(i, 0));
    double sum_sqs = square(x_ref(i, 0));
    for (int j = 1; j < i; ++j) {
      z(k++) = corr_free(x_ref(i, j) / std::sqrt(1.0 - sum_sqs));
      sum_sqs += square(x_ref(i, j));
    }
  }
  return z;
}

}  // namespace math
}  // namespace stan
#endif
