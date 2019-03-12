#ifndef STAN_MATH_PRIM_MAT_FUN_OPENCL_SYMMETRIC_EIGENSOLVER_HPP
#define STAN_MATH_PRIM_MAT_FUN_OPENCL_SYMMETRIC_EIGENSOLVER_HPP

#ifdef STAN_OPENCL

#include <queue>

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/multiply.hpp>
#include <stan/math/opencl/subtract.hpp>
#include <stan/math/opencl/add.hpp>
#include <stan/math/opencl/transpose.hpp>

#include <stan/math/opencl/tridiagonalization.hpp>
#include <stan/math/opencl/mrrr.hpp>

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

/**
 * Calculates eigenvalues and eigenvectors of a symmetric matrix.
 * @param A The matrix
 * @param eigenvals[out] Eigenvalues.
 * @param eigenvecs[out] Eigenvectors - one per column.
 */
void symmetric_eigensolver_cl(const Eigen::MatrixXd& A, Eigen::VectorXd& eigenvals, Eigen::MatrixXd& eigenvecs) {
  Eigen::MatrixXd packed;
  block_householder_tridiag_cl(A, packed);

  Eigen::VectorXd diag = packed.diagonal();
  Eigen::VectorXd subdiag = packed.diagonal(1);
  tridiagonal_eigensolver_cl(diag, subdiag, eigenvals, eigenvecs);

  block_apply_packed_Q_cl(packed, eigenvecs);
}

}
}
#endif
#endif
