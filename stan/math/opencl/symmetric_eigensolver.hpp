#ifndef STAN_MATH_PRIM_MAT_FUN_OPENCL_SYMMETRIC_EIGENSOLVER_HPP
#define STAN_MATH_PRIM_MAT_FUN_OPENCL_SYMMETRIC_EIGENSOLVER_HPP

#ifdef STAN_OPENCL

#include <Eigen/QR>
#include <iostream>
#include <queue>

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/multiply.hpp>
#include <stan/math/opencl/subtract.hpp>
#include <stan/math/opencl/add.hpp>
#include <stan/math/opencl/transpose.hpp>

#include <stan/math/opencl/tridiagonalization.hpp>
#include <stan/math/opencl/mrrr.hpp>

using namespace std;

#define TIME_IT
//#define SKIP_Q

namespace stan {
namespace math {

/**
 * Calculates eigenvalues and eigenvectors of a symmetric matrix.
 * @param A The matrix
 * @param eigenvals[out] Eigenvalues.b
 * @param eigenvecs[out] Eigenvectors.
 */
void symmetric_eigensolver(const Eigen::MatrixXd& A, Eigen::VectorXd& eigenvals, Eigen::MatrixXd& eigenvecs) {
  Eigen::MatrixXd packed;
  auto start = std::chrono::steady_clock::now();
  block_householder_tridiag4(A, packed);

  cout << "tridiag: "
       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
       << "ms" << endl;
  start = std::chrono::steady_clock::now();
  Eigen::VectorXd diag = packed.diagonal();
  Eigen::VectorXd subdiag = packed.diagonal(1);
  cout << "extract: "
       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
       << "ms" << endl;
  start = std::chrono::steady_clock::now();
  tridiagonal_eigensolver(diag, subdiag, eigenvals, eigenvecs);

  cout << "mrrr: "
       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
       << "ms" << endl;
  start = std::chrono::steady_clock::now();
  block_apply_packed_Q3(packed, eigenvecs);
  cout << "apply q: "
       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
       << "ms" << endl;
}

/**
 * Calculates eigenvalues and eigenvectors of a symmetric matrix.
 * @param A The matrix
 * @param eigenvals[out] Eigenvalues.b
 * @param eigenvecs[out] Eigenvectors.
 */
void symmetric_eigensolver_cl(const Eigen::MatrixXd& A, Eigen::VectorXd& eigenvals, Eigen::MatrixXd& eigenvecs) {
  Eigen::MatrixXd packed;
  auto start = std::chrono::steady_clock::now();
  block_householder_tridiag_cl2(A, packed);

  cout << "tridiag: "
       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
       << "ms" << endl;
  start = std::chrono::steady_clock::now();
  Eigen::VectorXd diag = packed.diagonal();
  Eigen::VectorXd subdiag = packed.diagonal(1);
  cout << "extract: "
       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
       << "ms" << endl;
  start = std::chrono::steady_clock::now();
  tridiagonal_eigensolver_cl(diag, subdiag, eigenvals, eigenvecs);

  cout << "mrrr: "
       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
       << "ms" << endl;
  start = std::chrono::steady_clock::now();
  block_apply_packed_Q_cl2(packed, eigenvecs);
  cout << "apply q: "
       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
       << "ms" << endl;
}

}
}
#endif //STAN_OPENCL
#endif //STAN_MATH_PRIM_MAT_FUN_OPENCL_EIGENDECOMPOSITION_HPP
