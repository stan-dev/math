#ifndef STAN_MATH_PRIM_MAT_FUN_SYMMETRIC_EIGENSOLVER_HPP
#define STAN_MATH_PRIM_MAT_FUN_SYMMETRIC_EIGENSOLVER_HPP

#include <stan/math/prim/mat/fun/tridiagonalization.hpp>
#include <stan/math/prim/mat/fun/mrrr.hpp>

#include <Eigen/Dense>
#include <queue>

namespace stan {
namespace math {

/**
 * Calculates eigenvalues and eigenvectors of a symmetric matrix.
 * @param A The matrix
 * @param[out] eigenvalues Eigenvalues.
 * @param[out] eigenvectors Eigenvectors - one per column.
 */
void symmetric_eigensolver(const Eigen::MatrixXd& A,
                           Eigen::VectorXd& eigenvalues,
                           Eigen::MatrixXd& eigenvectors) {
  Eigen::MatrixXd packed;
  internal::block_householder_tridiag(A, packed);
  Eigen::VectorXd diagonal = packed.diagonal();
  Eigen::VectorXd subdiagonal = packed.diagonal(1);
  internal::tridiagonal_eigensolver(diagonal, subdiagonal, eigenvalues,
                                    eigenvectors);
  internal::block_apply_packed_Q(packed, eigenvectors);
}

}  // namespace math
}  // namespace stan
#endif
