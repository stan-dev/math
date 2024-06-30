#ifndef STAN_MATH_PRIM_CONSTRAINT_CHOLESKY_FACTOR_FREE_HPP
#define STAN_MATH_PRIM_CONSTRAINT_CHOLESKY_FACTOR_FREE_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <cmath>
#include <stdexcept>

namespace stan {
namespace math {

/**
 * Return the unconstrained vector of parameters corresponding to
 * the specified Cholesky factor.  A Cholesky factor must be lower
 * triangular and have positive diagonal elements.
 *
 * @tparam T type of the Cholesky factor (must be derived from \c
 * Eigen::MatrixBase)
 * @param y Cholesky factor.
 * @return Unconstrained parameters for Cholesky factor.
 * @throw std::domain_error If the matrix is not a Cholesky factor.
 */
template <typename T, require_eigen_t<T>* = nullptr>
Eigen::Matrix<value_type_t<T>, Eigen::Dynamic, 1> cholesky_factor_free(
    const T& y) {
  using std::log;

  const auto& y_ref = to_ref(y);
  check_cholesky_factor("cholesky_factor_free", "y", y_ref);
  int M = y.rows();
  int N = y.cols();
  Eigen::Matrix<value_type_t<T>, Eigen::Dynamic, 1> x((N * (N + 1)) / 2
                                                      + (M - N) * N);
  int pos = 0;

  for (int m = 0; m < N; ++m) {
    x.segment(pos, m) = y_ref.row(m).head(m);
    pos += m;
    x.coeffRef(pos++) = log(y_ref.coeff(m, m));
  }

  for (int m = N; m < M; ++m) {
    x.segment(pos, N) = y_ref.row(m);
    pos += N;
  }
  return x;
}

/**
 * Overload of `cholesky_factor_free()` to untransform each matrix
 * in a standard vector.
 * @tparam T A standard vector with with a `value_type` which inherits from
 *  `Eigen::MatrixBase`.
 * @param x The standard vector to untransform.
 */
template <typename T, require_std_vector_t<T>* = nullptr>
auto cholesky_factor_free(const T& x) {
  return apply_vector_unary<T>::apply(
      x, [](auto&& v) { return cholesky_factor_free(v); });
}

}  // namespace math
}  // namespace stan

#endif
