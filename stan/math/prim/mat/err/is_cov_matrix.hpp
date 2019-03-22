#ifndef STAN_MATH_PRIM_MAT_ERR_IS_COV_MATRIX_HPP
#define STAN_MATH_PRIM_MAT_ERR_IS_COV_MATRIX_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/err/is_pos_definite.hpp>

namespace stan {
namespace math {
/**
 * Return <code>true</code> if the matrix is square or if the matrix
 * is 0x0, if the matrix is symmetric, if the matrix is positive
 * definite, or if no element of the matrix is <code>NaN</code>.
 * A valid covariance matrix is a square, symmetric matrix that is
 * positive definite.
 * @tparam T_y Type of scalar
 * @param y Matrix to test
 * @return <code>true</code> if the matrix is square or if the matrix
 *   is 0x0, if the matrix is symmetric, if the matrix is positive
 *   definite, or if no element of the matrix is <code>NaN</code>
 */
template <typename T_y>
inline bool is_cov_matrix(
    const Eigen::Matrix<T_y, Eigen::Dynamic, Eigen::Dynamic>& y) {
  return is_pos_definite(y);
}

}  // namespace math
}  // namespace stan
#endif
