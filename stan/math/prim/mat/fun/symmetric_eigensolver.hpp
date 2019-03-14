#ifndef STAN_MATH_PRIM_MAT_FUN_SYMMETRIC_EIGENSOLVER_HPP
#define STAN_MATH_PRIM_MAT_FUN_SYMMETRIC_EIGENSOLVER_HPP

#include <queue>

#include <Eigen/Dense>

#include <stan/math/prim/mat/fun/tridiagonalization.hpp>
#include <stan/math/prim/mat/fun/mrrr.hpp>

#include <stan/math/opencl/copy.hpp>
#include <stan/math/opencl/tridiagonalization.hpp>
#include <stan/math/opencl/mrrr.hpp>

namespace stan {
namespace math {

/**
 * Calculates eigenvalues and eigenvectors of a symmetric matrix.
 * @param A The matrix
 * @param eigenvalues[out] Eigenvalues.
 * @param eigenvectors[out] Eigenvectors - one per column.
 */
void symmetric_eigensolver(const Eigen::MatrixXd& A, Eigen::VectorXd& eigenvalues, Eigen::MatrixXd& eigenvectors) {
  Eigen::MatrixXd packed;
#ifdef STAN_OPENCL
  internal::block_householder_tridiag_cl(A, packed);
#else
  internal::block_householder_tridiag(A, packed);
#endif
  Eigen::VectorXd diagonal = packed.diagonal();
  Eigen::VectorXd subdiagonal = packed.diagonal(1);
#ifdef STAN_OPENCL
  internal::tridiagonal_eigensolver_cl(diagonal, subdiagonal, eigenvalues, eigenvectors);
  internal::block_apply_packed_Q_cl(packed, eigenvectors);
#else
  internal::tridiagonal_eigensolver(diagonal, subdiagonal, eigenvalues, eigenvectors);
  internal::block_apply_packed_Q(packed, eigenvectors);
#endif
}

}
}
#endif
