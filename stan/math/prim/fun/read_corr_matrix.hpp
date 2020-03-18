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
 * @tparam T_CPCs type of the array (must be derived from \c Eigen::ArrayBase
 * and have one compile-time dimension equal to 1)
 * @param CPCs The (K choose 2) canonical partial correlations in (-1, 1).
 * @param K Dimensionality of correlation matrix.
 * @return Cholesky factor of correlation matrix for specified
 * canonical partial correlations.
 */
template <typename T_CPCs, require_eigen_vector_t<T_CPCs>* = nullptr>
Eigen::Matrix<value_type_t<T_CPCs>, Eigen::Dynamic, Eigen::Dynamic>
read_corr_matrix(const T_CPCs& CPCs, size_t K) {
  if (K == 0) {
    return {};
  }

  return multiply_lower_tri_self_transpose(read_corr_L(CPCs, K));
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
 * @tparam T_CPCs type of the array (must be derived from \c Eigen::ArrayBase
 * and have one compile-time dimension equal to 1)
 * @param CPCs The (K choose 2) canonical partial correlations in
 * (-1, 1).
 * @param K Dimensionality of correlation matrix.
 * @param log_prob Reference to variable to increment with the log
 * Jacobian determinant.
 * @return Correlation matrix for specified partial correlations.
 */
template <typename T_CPCs, require_eigen_vector_t<T_CPCs>* = nullptr>
Eigen::Matrix<value_type_t<T_CPCs>, Eigen::Dynamic, Eigen::Dynamic>
read_corr_matrix(const T_CPCs& CPCs, size_t K, value_type_t<T_CPCs>& log_prob) {
  if (K == 0) {
    return {};
  }

  return multiply_lower_tri_self_transpose(read_corr_L(CPCs, K, log_prob));
}

}  // namespace math
}  // namespace stan

#endif
