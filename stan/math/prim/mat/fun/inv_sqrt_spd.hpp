#ifndef STAN_MATH_PRIM_MAT_FUN_INV_SQRT_SPD_HPP
#define STAN_MATH_PRIM_MAT_FUN_INV_SQRT_SPD_HPP

#include <stan/math/prim/arr/err/check_nonzero_size.hpp>
#include <stan/math/prim/mat/err/check_symmetric.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Returns the inverse symmetric square root of the specified
 * symmetric and positive definite matrix, which is defined as
 * V * diag_matrix(inv_sqrt(d)) * V'
 * where V is a matrix of eigenvectors of the input matrix and
 * d is a vector of eigenvalue of the input matrix. If the
 * input matrix is not actually positive definite, then the
 * output matrix will have NaNs in it.
 *
 * @param[in] m Specified symmetric positive definite matrix.
 * @return Inverse symmetric square root of the matrix.
 * @throws If matrix is of size zero or asymmetric
 */

template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> inv_sqrt_spd(
    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &m) {
  check_nonzero_size("inv_sqrt_spd", "m", m);
  check_symmetric("inv_sqrt_spd", "m", m);

  Eigen::SelfAdjointEigenSolver<
      Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>
      solver(m);
  return solver.operatorInverseSqrt();
}

}  // namespace math
}  // namespace stan
#endif
