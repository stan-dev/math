#ifndef STAN_MATH_PRIM_FUN_CHOLESKY_FACTOR_CONSTRAIN_HPP
#define STAN_MATH_PRIM_FUN_CHOLESKY_FACTOR_CONSTRAIN_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace stan {
namespace math {

/**
 * Return the Cholesky factor of the specified size read from the
 * specified vector.  A total of (N choose 2) + N + (M - N) * N
 * elements are required to read an M by N Cholesky factor.
 *
 * @tparam T type of the vector
 * @param x Vector of unconstrained values
 * @return Cholesky factor
 */
template <typename T, require_eigen_vector_t<T>* = nullptr>
Eigen::Matrix<value_type_t<T>, Eigen::Dynamic, Eigen::Dynamic> cholesky_factor_constrain(
    const T& x, int M, int N) {
  using std::exp;
  using T_scalar = value_type_t<T>;
  check_greater_or_equal("cholesky_factor_constrain",
                         "num rows (must be greater or equal to num cols)", M,
                         N);
  check_size_match("cholesky_factor_constrain", "x.size()", x.size(),
                   "((N * (N + 1)) / 2 + (M - N) * N)",
                   ((N * (N + 1)) / 2 + (M - N) * N));
  Eigen::Matrix<T_scalar, Eigen::Dynamic, Eigen::Dynamic> y(M, N);
  T_scalar zero(0);
  int pos = 0;

  const Eigen::Ref<const plain_type_t<T>>& x_ref = x;
  for (int m = 0; m < N; ++m) {
    for (int n = 0; n < m; ++n) {
      y(m, n) = x_ref(pos++);
    }
    y(m, m) = exp(x_ref(pos++));
    for (int n = m + 1; n < N; ++n) {
      y(m, n) = zero;
    }
  }

  for (int m = N; m < M; ++m) {
    for (int n = 0; n < N; ++n) {
      y(m, n) = x_ref(pos++);
    }
  }
  return y;
}

/**
 * Return the Cholesky factor of the specified size read from the
 * specified vector and increment the specified log probability
 * reference with the log Jacobian adjustment of the transform.  A total
 * of (N choose 2) + N + N * (M - N) free parameters are required to read
 * an M by N Cholesky factor.
 *
 * @tparam T type of elements in the matrix
 * @param x Vector of unconstrained values
 * @param lp Log probability that is incremented with the log Jacobian
 * @return Cholesky factor
 */
template <typename T, require_eigen_vector_t<T>* = nullptr>
Eigen::Matrix<value_type_t<T>, Eigen::Dynamic, Eigen::Dynamic> cholesky_factor_constrain(
    const T& x, int M, int N, value_type_t<T>& lp) {
  check_size_match("cholesky_factor_constrain", "x.size()", x.size(),
                   "((N * (N + 1)) / 2 + (M - N) * N)",
                   ((N * (N + 1)) / 2 + (M - N) * N));
  int pos = 0;
  std::vector<value_type_t<T>> log_jacobians(N);
  const Eigen::Ref<const plain_type_t<T>>& x_ref = x;
  for (int n = 0; n < N; ++n) {
    pos += n;
    log_jacobians[n] = x_ref(pos++);
  }
  lp += sum(log_jacobians);
  return cholesky_factor_constrain(x_ref, M, N);
}

}  // namespace math
}  // namespace stan

#endif
