#ifndef STAN_MATH_PRIM_CONSTRAINT_CHOLESKY_CORR_FREE_HPP
#define STAN_MATH_PRIM_CONSTRAINT_CHOLESKY_CORR_FREE_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/constraint/corr_free.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T, require_eigen_t<T>* = nullptr>
auto cholesky_corr_free(const T& x) {
  using Eigen::Dynamic;
  using Eigen::Matrix;

  check_square("cholesky_corr_free", "x", x);
  // should validate lower-triangular, unit lengths

  const auto& x_ref = to_ref(x);
  int K = (x.rows() * (x.rows() - 1)) / 2;
  Matrix<value_type_t<T>, Dynamic, 1> z(K);
  int k = 0;
  for (int i = 1; i < x.rows(); ++i) {
    z.coeffRef(k++) = corr_free(x_ref.coeff(i, 0));
    double sum_sqs = square(x_ref.coeff(i, 0));
    for (int j = 1; j < i; ++j) {
      z.coeffRef(k++) = corr_free(x_ref.coeff(i, j) / std::sqrt(1.0 - sum_sqs));
      sum_sqs += square(x_ref.coeff(i, j));
    }
  }
  return z;
}

/**
 * Overload of `cholesky_corr_free()` to untransform each matrix
 * in a standard vector.
 * @tparam T A standard vector with with a `value_type` which inherits from
 *  `Eigen::MatrixBase`.
 * @param x The standard vector to untransform.
 */
template <typename T, require_std_vector_t<T>* = nullptr>
auto cholesky_corr_free(const T& x) {
  return apply_vector_unary<T>::apply(
      x, [](auto&& v) { return cholesky_corr_free(v); });
}

}  // namespace math
}  // namespace stan
#endif
