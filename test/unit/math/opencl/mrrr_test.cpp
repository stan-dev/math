#ifdef STAN_OPENCL
#include <iostream>
#include <stan/math/opencl/mrrr.hpp>
#include <test/unit/math/expect_near_rel.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(MathMatrix, tridiag_eigensolver_trivial) {
  Eigen::VectorXd diag(3), subdiag(2), eigenvals;
  diag << 1.5, 1.2, -2;
  subdiag << 0, 0;
  stan::math::matrix_cl<double> diag_cl(diag);
  stan::math::matrix_cl<double> subdiag_cl(subdiag);
  stan::math::matrix_cl<double> eigenvals_cl;
  stan::math::matrix_cl<double> eigenvecs_cl;

  stan::math::internal::tridiagonal_eigensolver_cl(diag_cl, subdiag_cl,
                                                   eigenvals_cl, eigenvecs_cl);
  EXPECT_NEAR_REL(stan::math::from_matrix_cl(eigenvals_cl), diag);
  EXPECT_NEAR_REL(stan::math::from_matrix_cl(eigenvecs_cl),
                  Eigen::MatrixXd::Identity(3, 3));
}

TEST(MathMatrix, tridiag_eigensolver_small) {
  int size = 7;
  Eigen::VectorXd diag = Eigen::VectorXd::Constant(size, 1);
  Eigen::VectorXd subdiag = Eigen::VectorXd::Constant(size - 1, 0.5);
  stan::math::matrix_cl<double> diag_cl(diag);
  stan::math::matrix_cl<double> subdiag_cl(subdiag);
  stan::math::matrix_cl<double> eigenvals_cl;
  stan::math::matrix_cl<double> eigenvecs_cl;

  stan::math::internal::tridiagonal_eigensolver_cl(diag_cl, subdiag_cl,
                                                   eigenvals_cl, eigenvecs_cl);

  Eigen::VectorXd eigenvals
      = stan::math::from_matrix_cl<Eigen::VectorXd>(eigenvals_cl);
  Eigen::MatrixXd eigenvecs = stan::math::from_matrix_cl(eigenvecs_cl);
  Eigen::MatrixXd t = Eigen::MatrixXd::Constant(size, size, 0);
  t.diagonal() = diag;
  t.diagonal(1) = subdiag;
  t.diagonal(-1) = subdiag;
  EXPECT_NEAR_REL(diag.sum(), eigenvals.sum());
  EXPECT_MATRIX_NEAR(eigenvecs * eigenvecs.transpose(),
                     Eigen::MatrixXd::Identity(size, size), 1e-12);
  EXPECT_MATRIX_NEAR(t * eigenvecs, eigenvecs * eigenvals.asDiagonal(), 1e-12);
}

TEST(MathMatrix, tridiag_eigensolver_large) {
  int size = 2000;
  Eigen::VectorXd diag = Eigen::VectorXd::Random(size);
  Eigen::VectorXd subdiag = Eigen::VectorXd::Random(size - 1);
  subdiag[12] = 0;
  subdiag[120] = 0;
  subdiag[121] = 0;

  stan::math::matrix_cl<double> diag_cl(diag);
  stan::math::matrix_cl<double> subdiag_cl(subdiag);
  stan::math::matrix_cl<double> eigenvals_cl;
  stan::math::matrix_cl<double> eigenvecs_cl;

  stan::math::internal::tridiagonal_eigensolver_cl(diag_cl, subdiag_cl,
                                                   eigenvals_cl, eigenvecs_cl);

  Eigen::VectorXd eigenvals
      = stan::math::from_matrix_cl<Eigen::VectorXd>(eigenvals_cl);
  Eigen::MatrixXd eigenvecs = stan::math::from_matrix_cl(eigenvecs_cl);

  Eigen::MatrixXd t = Eigen::MatrixXd::Constant(size, size, 0);
  t.diagonal() = diag;
  t.diagonal(1) = subdiag;
  t.diagonal(-1) = subdiag;
  EXPECT_NEAR_REL(diag.sum(), eigenvals.sum());
  EXPECT_MATRIX_NEAR(eigenvecs * eigenvecs.transpose(),
                     Eigen::MatrixXd::Identity(size, size), 1e-12);
  EXPECT_MATRIX_NEAR(t * eigenvecs, eigenvecs * eigenvals.asDiagonal(), 1e-12);
}

TEST(MathMatrix, tridiag_eigensolver_large_fixed) {
  int size = 2000;
  Eigen::VectorXd diag = Eigen::VectorXd::Constant(size, 0.4);
  Eigen::VectorXd subdiag = Eigen::VectorXd::Constant(size - 1, 0.4);

  stan::math::matrix_cl<double> diag_cl(diag);
  stan::math::matrix_cl<double> subdiag_cl(subdiag);
  stan::math::matrix_cl<double> eigenvals_cl;
  stan::math::matrix_cl<double> eigenvecs_cl;

  stan::math::internal::tridiagonal_eigensolver_cl(diag_cl, subdiag_cl,
                                                   eigenvals_cl, eigenvecs_cl);

  Eigen::VectorXd eigenvals
      = stan::math::from_matrix_cl<Eigen::VectorXd>(eigenvals_cl);
  Eigen::MatrixXd eigenvecs = stan::math::from_matrix_cl(eigenvecs_cl);

  Eigen::MatrixXd t = Eigen::MatrixXd::Constant(size, size, 0);
  t.diagonal() = diag;
  t.diagonal(1) = subdiag;
  t.diagonal(-1) = subdiag;
  EXPECT_NEAR_REL(diag.sum(), eigenvals.sum());
  EXPECT_MATRIX_NEAR(eigenvecs * eigenvecs.transpose(),
                     Eigen::MatrixXd::Identity(size, size), 1e-12);
  EXPECT_MATRIX_NEAR(t * eigenvecs, eigenvecs * eigenvals.asDiagonal(), 1e-12);
}

// COMMENTED OUT DUE TO INFINITE LOOP BEHAVIOR WITH CUDA 12
// /**
//  * Generates a glued Wilkinson matrix - a known hard case for MRRR.
//  * @param m each Wilkinson matrix is of size `2*m+1`
//  * @param n number of Wilkinson matrices to use
//  * @param glue amount og glue
//  * @param[out] diag diagonal of the resulting matrix
//  * @param[out] subdiag subdiagonal of the resulting matrix
//  */
// void get_glued_wilkinson(int m, int n, double glue, Eigen::VectorXd& diag,
//                          Eigen::VectorXd& subdiag) {
//   int w_size = 2 * m + 1;
//   int size = w_size * n;
//   diag.resize(size);
//   subdiag.resize(size - 1);
//   for (int start = 0; start < size; start += w_size) {
//     for (int i = 0; i < m; i++) {
//       int value = m - i;
//       diag[start + i] = value;
//       diag[start + w_size - i - 1] = value;
//     }
//     diag[start + m] = 0;
//     subdiag.segment(start, w_size - 1)
//         = Eigen::VectorXd::Constant(w_size - 1, 1);
//     if (start > 0) {
//       subdiag[start - 1] = glue;
//       diag[start] += glue;
//       diag[start - 1] += glue;
//     }
//   }
// }

// TEST(MathMatrix, tridiag_eigensolver_large_wilkinson) {
//   int n = 5;
//   int m = 100;
//   int size = (2 * m + 1) * n;
//   Eigen::VectorXd diag;
//   Eigen::VectorXd subdiag;
//   get_glued_wilkinson(m, n, 1e-13, diag, subdiag);

//   stan::math::matrix_cl<double> diag_cl(diag);
//   stan::math::matrix_cl<double> subdiag_cl(subdiag);
//   stan::math::matrix_cl<double> eigenvals_cl;
//   stan::math::matrix_cl<double> eigenvecs_cl;

//   stan::math::internal::tridiagonal_eigensolver_cl(diag_cl, subdiag_cl,
//                                                    eigenvals_cl,
//                                                    eigenvecs_cl);

//   Eigen::VectorXd eigenvals
//       = stan::math::from_matrix_cl<Eigen::VectorXd>(eigenvals_cl);
//   Eigen::MatrixXd eigenvecs = stan::math::from_matrix_cl(eigenvecs_cl);

//   Eigen::MatrixXd t = Eigen::MatrixXd::Constant(size, size, 0);
//   t.diagonal() = diag;
//   t.diagonal(1) = subdiag;
//   t.diagonal(-1) = subdiag;
//   EXPECT_NEAR_REL(diag.sum(), eigenvals.sum());
//   EXPECT_MATRIX_NEAR(eigenvecs * eigenvecs.transpose(),
//                      Eigen::MatrixXd::Identity(size, size), 1e-11);
//   EXPECT_MATRIX_NEAR(t * eigenvecs, eigenvecs * eigenvals.asDiagonal(),
//   1e-12);
// }

#endif
