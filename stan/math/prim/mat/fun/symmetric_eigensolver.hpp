#ifndef STAN_MATH_PRIM_MAT_FUN_SYMMETRIC_EIGENSOLVER_HPP
#define STAN_MATH_PRIM_MAT_FUN_SYMMETRIC_EIGENSOLVER_HPP

#include <queue>

#include <Eigen/Dense>

#include <stan/math/prim/mat/fun/tridiagonalization.hpp>
#include <stan/math/prim/mat/fun/mrrr.hpp>

namespace stan {
namespace math {

/**
 * Calculates eigenvalues and eigenvectors of a symmetric matrix.
 * @param A The matrix
 * @param eigenvals[out] Eigenvalues.
 * @param eigenvecs[out] Eigenvectors - one per column.
 */
void symmetric_eigensolver(const Eigen::MatrixXd& A, Eigen::VectorXd& eigenvals, Eigen::MatrixXd& eigenvecs) {
  Eigen::MatrixXd packed;
  block_householder_tridiag(A, packed);
  Eigen::VectorXd diag = packed.diagonal();
  Eigen::VectorXd subdiag = packed.diagonal(1);
  tridiagonal_eigensolver(diag, subdiag, eigenvals, eigenvecs);
  block_apply_packed_Q(packed, eigenvecs);
}

}
}
#endif
