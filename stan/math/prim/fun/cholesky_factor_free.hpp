#ifndef STAN_MATH_PRIM_FUN_CHOLESKY_FACTOR_FREE_HPP
#define STAN_MATH_PRIM_FUN_CHOLESKY_FACTOR_FREE_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <cmath>
#include <stdexcept>

namespace stan {
namespace math {

/**
 * Return the unconstrained vector of parameters correspdonding to
 * the specified Cholesky factor.  A Cholesky factor must be lower
 * triangular and have positive diagonal elements.
 *
 * @tparam T type of elements in the Cholesky factor
 * @param y Cholesky factor.
 * @return Unconstrained parameters for Cholesky factor.
 * @throw std::domain_error If the matrix is not a Cholesky factor.
 */
template <typename T, require_eigen_t<T>* = nullptr>
Eigen::Matrix<value_type_t<T>, Eigen::Dynamic, 1> cholesky_factor_free(
    const T& y) {
  using std::log;

  const Eigen::Ref<const plain_type_t<T>>& y_ref = y;
  check_cholesky_factor("cholesky_factor_free", "y", y_ref);
  int M = y.rows();
  int N = y.cols();
  Eigen::Matrix<value_type_t<T>, Eigen::Dynamic, 1> x((N * (N + 1)) / 2
                                                      + (M - N) * N);
  int pos = 0;

  for (int m = 0; m < N; ++m) {
    for (int n = 0; n < m; ++n) {
      x(pos++) = y_ref(m, n);
    }
    x(pos++) = log(y_ref(m, m));
  }

  for (int m = N; m < M; ++m) {
    for (int n = 0; n < N; ++n) {
      x(pos++) = y_ref(m, n);
    }
  }
  return x;
}

}  // namespace math
}  // namespace stan

#endif
