#ifndef STAN_MATH_PRIM_FUN_READ_CORR_MATRIX_HPP
#define STAN_MATH_PRIM_FUN_READ_CORR_MATRIX_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/read_corr_L.hpp>
#include <stan/math/prim/fun/multiply_lower_tri_self_transpose.hpp>

namespace stan {
namespace math {

/**
 * Return the correlation matrix of the specified dimensionality
 * corresponding to the specified canonical partial correlations.
 *
 * <p>See <code>read_corr_matrix(Array, size_t, T)</code>
 * for more information.
 *
 * @tparam T type of elements in the array
 * @param CPCs The (K choose 2) canonical partial correlations in (-1, 1).
 * @param K Dimensionality of correlation matrix.
 * @return Cholesky factor of correlation matrix for specified
 * canonical partial correlations.
 */
template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> read_corr_matrix(
    const Eigen::Array<T, Eigen::Dynamic, 1>& CPCs, size_t K) {
  if (K == 0) {
    return {};
  }

  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> L = read_corr_L(CPCs, K);
  return multiply_lower_tri_self_transpose(L);
}

/**
 * Return the correlation matrix of the specified dimensionality
 * corresponding to the specified canonical partial correlations,
 * incrementing the specified scalar reference with the log
 * absolute determinant of the Jacobian of the transformation.
 *
 * It is usually preferable to utilize the version that returns
 * the Cholesky factor of the correlation matrix rather than the
 * correlation matrix itself in statistical calculations.
 *
 * @tparam T type of elements in the array
 * @param CPCs The (K choose 2) canonical partial correlations in
 * (-1, 1).
 * @param K Dimensionality of correlation matrix.
 * @param log_prob Reference to variable to increment with the log
 * Jacobian determinant.
 * @return Correlation matrix for specified partial correlations.
 */
template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> read_corr_matrix(
    const Eigen::Array<T, Eigen::Dynamic, 1>& CPCs, size_t K, T& log_prob) {
  if (K == 0) {
    return {};
  }

  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> L
      = read_corr_L(CPCs, K, log_prob);
  return multiply_lower_tri_self_transpose(L);
}

}  // namespace math
}  // namespace stan

#endif
