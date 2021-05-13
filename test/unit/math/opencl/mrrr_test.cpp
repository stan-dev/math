
#include <iostream>
#include <stan/math/opencl/mrrr.hpp>
#include <gtest/gtest.h>

TEST(MathMatrix, tridiag_eigensolver_trivial) {
  Eigen::VectorXd diag(3), subdiag(2), eigenvals;
  diag << 1.5, 1.2, -2;
  subdiag << 0, 0;

  Eigen::MatrixXd eigenvecs(3, 3);
  stan::math::internal::tridiagonal_eigensolver_cl(diag, subdiag, eigenvals,
                                                   eigenvecs);
  EXPECT_TRUE(eigenvals.isApprox(diag));
  EXPECT_TRUE(eigenvecs.isApprox(Eigen::MatrixXd::Identity(3, 3)));
}

TEST(MathMatrix, tridiag_eigensolver_small) {
  int size = 7;
  srand(0);  // ensure test repeatability
  Eigen::VectorXd diag = Eigen::VectorXd::Random(size);
  Eigen::VectorXd subdiag = Eigen::VectorXd::Random(size - 1);
  subdiag[2] = 0;

  Eigen::VectorXd eigenvals;
  Eigen::MatrixXd eigenvecs;
  stan::math::internal::tridiagonal_eigensolver_cl(diag, subdiag, eigenvals,
                                                   eigenvecs);

  Eigen::MatrixXd t = Eigen::MatrixXd::Constant(size, size, 0);
  t.diagonal() = diag;
  t.diagonal(1) = subdiag;
  t.diagonal(-1) = subdiag;
  EXPECT_NEAR(diag.sum(), eigenvals.sum(), 1e-11);
  EXPECT_TRUE((eigenvecs * eigenvecs.transpose())
                  .isApprox(Eigen::MatrixXd::Identity(size, size)));
  EXPECT_TRUE((t * eigenvecs).isApprox(eigenvecs * eigenvals.asDiagonal()));
}

TEST(MathMatrix, tridiag_eigensolver_large) {
  int size = 2000;
  srand(time(0));  // ensure test repeatability
  Eigen::VectorXd diag = Eigen::VectorXd::Random(size);
  Eigen::VectorXd subdiag = Eigen::VectorXd::Random(size - 1);
  subdiag[12] = 0;
  subdiag[120] = 0;
  subdiag[121] = 0;
  //  Eigen::VectorXd diag = Eigen::VectorXd::Constant(size, 0.4);
  //  Eigen::VectorXd subdiag = Eigen::VectorXd::Constant(size - 1, 0.4);

  Eigen::VectorXd eigenvals;
  Eigen::MatrixXd eigenvecs;
  //  stan::math::internal::tridiagonal_eigensolver_cl(diag, subdiag, eigenvals,
  //                                                   eigenvecs);
  stan::math::internal::tridiagonal_eigensolver(diag, subdiag, eigenvals,
                                                eigenvecs);
  //  Eigen::EigenSolver

  Eigen::MatrixXd t = Eigen::MatrixXd::Constant(size, size, 0);
  t.diagonal() = diag;
  t.diagonal(1) = subdiag;
  t.diagonal(-1) = subdiag;
  EXPECT_NEAR(diag.sum(), eigenvals.sum(), 1e-9);
  EXPECT_TRUE((eigenvecs * eigenvecs.transpose())
                  .isApprox(Eigen::MatrixXd::Identity(size, size), 1e-9));
  EXPECT_TRUE(
      (t * eigenvecs).isApprox(eigenvecs * eigenvals.asDiagonal(), 1e-9));
  //  std::cout << diag << std::endl << std::endl;
  //  std::cout << subdiag << std::endl << std::endl;
  //  std::cout << eigenvecs << std::endl << std::endl;
  //  std::cout << eigenvals << std::endl << std::endl;
  //  std::cout << eigenvecs * eigenvecs.transpose() << std::endl << std::endl;
  //  std::cout << t * eigenvecs  << std::endl << std::endl;
  //  std::cout << eigenvecs * eigenvals.asDiagonal() << std::endl << std::endl;
  //  std::cout << t * eigenvecs - eigenvecs * eigenvals.asDiagonal() <<
  //  std::endl << std::endl;
  std::cout << "trace rel err: " << abs(diag.sum() - eigenvals.sum())
            << std::endl;
  std::cout << "vec ortho err: "
            << ((eigenvecs * eigenvecs.transpose())
                - Eigen::MatrixXd::Identity(size, size))
                   .array()
                   .abs()
                   .maxCoeff()
            << std::endl;
  std::cout
      << "eigen eq err: "
      << ((t * eigenvecs - eigenvecs * eigenvals.asDiagonal()).array().abs())
             .array()
             .abs()
             .maxCoeff()
      << std::endl;
}

TEST(MathMatrix, tridiag_eigensolver_large_fixed) {
  int size = 2000;
  srand(time(0));  // ensure test repeatability
                   //  Eigen::VectorXd diag = Eigen::VectorXd::Random(size);
  //  Eigen::VectorXd subdiag = Eigen::VectorXd::Random(size - 1);
  //  subdiag[12] = 0;
  //  subdiag[120] = 0;
  //  subdiag[121] = 0;
  Eigen::VectorXd diag = Eigen::VectorXd::Constant(size, 0.4);
  Eigen::VectorXd subdiag = Eigen::VectorXd::Constant(size - 1, 0.4);

  Eigen::VectorXd eigenvals;
  Eigen::MatrixXd eigenvecs;
  //  stan::math::internal::tridiagonal_eigensolver_cl(diag, subdiag, eigenvals,
  //                                                   eigenvecs);
  stan::math::internal::tridiagonal_eigensolver(diag, subdiag, eigenvals,
                                                eigenvecs);
  //  Eigen::EigenSolver

  Eigen::MatrixXd t = Eigen::MatrixXd::Constant(size, size, 0);
  t.diagonal() = diag;
  t.diagonal(1) = subdiag;
  t.diagonal(-1) = subdiag;
  EXPECT_NEAR(diag.sum(), eigenvals.sum(), 1e-9);
  EXPECT_TRUE((eigenvecs * eigenvecs.transpose())
                  .isApprox(Eigen::MatrixXd::Identity(size, size), 1e-9));
  EXPECT_TRUE(
      (t * eigenvecs).isApprox(eigenvecs * eigenvals.asDiagonal(), 1e-9));
  //  std::cout << diag << std::endl << std::endl;
  //  std::cout << subdiag << std::endl << std::endl;
  //  std::cout << eigenvecs << std::endl << std::endl;
  //  std::cout << eigenvals << std::endl << std::endl;
  //  std::cout << eigenvecs * eigenvecs.transpose() << std::endl << std::endl;
  //  std::cout << t * eigenvecs  << std::endl << std::endl;
  //  std::cout << eigenvecs * eigenvals.asDiagonal() << std::endl << std::endl;
  //  std::cout << t * eigenvecs - eigenvecs * eigenvals.asDiagonal() <<
  //  std::endl << std::endl;
  std::cout << "trace rel err: " << abs(diag.sum() - eigenvals.sum())
            << std::endl;
  std::cout << "vec ortho err: "
            << ((eigenvecs * eigenvecs.transpose())
                - Eigen::MatrixXd::Identity(size, size))
                   .array()
                   .abs()
                   .maxCoeff()
            << std::endl;
  std::cout
      << "eigen eq err: "
      << ((t * eigenvecs - eigenvecs * eigenvals.asDiagonal()).array().abs())
             .array()
             .abs()
             .maxCoeff()
      << std::endl;
}

/**
 * Generates a glued Wilkinson matrix - a known hard case for MRRR.
 * @param m each Wilkinson matrix is of size `2*m+1`
 * @param n number of Wilkinson matrices to use
 * @param glue amount og glue
 * @param[out] diag diagonal of the resulting matrix
 * @param[out] subdiag subdiagonal of the resulting matrix
 */
void get_glued_wilkinson(int m, int n, double glue, Eigen::VectorXd& diag,
                         Eigen::VectorXd& subdiag) {
  int w_size = 2 * m + 1;
  int size = w_size * n;
  diag.resize(size);
  subdiag.resize(size - 1);
  for (int start = 0; start < size; start += w_size) {
    for (int i = 0; i < m; i++) {
      int value = m - i;
      diag[start + i] = value;
      diag[start + w_size - i - 1] = value;
    }
    diag[start + m] = 0;
    subdiag.segment(start, w_size - 1)
        = Eigen::VectorXd::Constant(w_size - 1, 1);
    if (start > 0) {
      subdiag[start - w_size] = glue;
      diag[start] += glue;
      diag[start - 1] += glue;
    }
  }
}

TEST(MathMatrix, tridiag_eigensolver_large_wilkinson) {
  int n = 5;
  int m = 100;
  int size = (2 * m + 1) * n;
  srand(time(0));  // ensure test repeatability
                   //  srand(0);  // ensure test repeatability
                   //  Eigen::VectorXd diag = Eigen::VectorXd::Random(size);
  //  Eigen::VectorXd subdiag = Eigen::VectorXd::Random(size - 1);
  //  subdiag[12] = 0;
  //  subdiag[120] = 0;
  //  subdiag[121] = 0;
  Eigen::VectorXd diag;
  Eigen::VectorXd subdiag;
  get_glued_wilkinson(m, n, 1e-15, diag, subdiag);

  Eigen::VectorXd eigenvals;
  Eigen::MatrixXd eigenvecs;
  //  stan::math::internal::tridiagonal_eigensolver_cl(diag, subdiag, eigenvals,
  //                                                   eigenvecs);
  stan::math::internal::tridiagonal_eigensolver(diag, subdiag, eigenvals,
                                                eigenvecs);
  //  Eigen::EigenSolver

  Eigen::MatrixXd t = Eigen::MatrixXd::Constant(size, size, 0);
  t.diagonal() = diag;
  t.diagonal(1) = subdiag;
  t.diagonal(-1) = subdiag;
  EXPECT_NEAR(diag.sum(), eigenvals.sum(), 1e-9);
  EXPECT_TRUE((eigenvecs * eigenvecs.transpose())
                  .isApprox(Eigen::MatrixXd::Identity(size, size), 1e-9));
  EXPECT_TRUE(
      (t * eigenvecs).isApprox(eigenvecs * eigenvals.asDiagonal(), 1e-9));
  //  std::cout << diag << std::endl << std::endl;
  //  std::cout << subdiag << std::endl << std::endl;
  //  std::cout << eigenvecs << std::endl << std::endl;
  //  std::cout << eigenvals << std::endl << std::endl;
  //  std::cout << eigenvecs * eigenvecs.transpose() << std::endl << std::endl;
  //  std::cout << t * eigenvecs  << std::endl << std::endl;
  //  std::cout << eigenvecs * eigenvals.asDiagonal() << std::endl << std::endl;
  //  std::cout << t * eigenvecs - eigenvecs * eigenvals.asDiagonal() <<
  //  std::endl << std::endl;
  std::cout << "trace rel err: "
            << abs(diag.sum() - eigenvals.sum()) / abs(diag.sum()) << std::endl;
  std::cout << "vec ortho err: "
            << ((eigenvecs * eigenvecs.transpose())
                - Eigen::MatrixXd::Identity(size, size))
                   .array()
                   .abs()
                   .maxCoeff()
            << std::endl;
  std::cout
      << "eigen eq rel err: "
      << (((t * eigenvecs - eigenvecs * eigenvals.asDiagonal()).array().abs())
              .array()
              .abs()
              .colwise()
          / eigenvals.array())
             .maxCoeff()
      << std::endl;
}
