#ifdef STAN_OPENCL
#include <iostream>
#include <stan/math/opencl/mrrr.hpp>
#include <test/unit/math/expect_near_rel.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(MathMatrix, symmetric_eigensolver_small) {
  int size = 7;
  Eigen::MatrixXd A = Eigen::MatrixXd::Random(size, size);
  A += A.transpose().eval();

  stan::math::matrix_cl<double> A_cl(A);
  stan::math::matrix_cl<double> eigenvals_cl;
  stan::math::matrix_cl<double> eigenvecs_cl;

  stan::math::symmetric_eigensolver(A_cl, eigenvals_cl, eigenvecs_cl);

  Eigen::VectorXd eigenvals
      = stan::math::from_matrix_cl<Eigen::VectorXd>(eigenvals_cl);
  Eigen::MatrixXd eigenvecs = stan::math::from_matrix_cl(eigenvecs_cl);
  EXPECT_NEAR_REL(A.diagonal().sum(), eigenvals.sum());
  EXPECT_MATRIX_NEAR(eigenvecs * eigenvecs.transpose(),
                     Eigen::MatrixXd::Identity(size, size), 1e-12);
  EXPECT_MATRIX_NEAR(A * eigenvecs, eigenvecs * eigenvals.asDiagonal(), 1e-12);
}

TEST(MathMatrix, symmetric_eigensolver_large) {
  int size = 2000;
  Eigen::MatrixXd A = Eigen::MatrixXd::Random(size, size);
  A += A.transpose().eval();

  stan::math::matrix_cl<double> A_cl(A);
  stan::math::matrix_cl<double> eigenvals_cl;
  stan::math::matrix_cl<double> eigenvecs_cl;

  stan::math::symmetric_eigensolver(A_cl, eigenvals_cl, eigenvecs_cl);

  Eigen::VectorXd eigenvals
      = stan::math::from_matrix_cl<Eigen::VectorXd>(eigenvals_cl);
  Eigen::MatrixXd eigenvecs = stan::math::from_matrix_cl(eigenvecs_cl);

  EXPECT_NEAR_REL(A.diagonal().sum(), eigenvals.sum());
  EXPECT_MATRIX_NEAR(eigenvecs * eigenvecs.transpose(),
                     Eigen::MatrixXd::Identity(size, size), 1e-12);
  EXPECT_MATRIX_NEAR(A * eigenvecs, eigenvecs * eigenvals.asDiagonal(), 1e-12);
}

TEST(MathMatrix, symmetric_eigensolver_prim_small) {
  for (int size = 1; size < 7; size++) {
    Eigen::MatrixXd A = Eigen::MatrixXd::Random(size, size);
    A += A.transpose().eval();

    stan::math::matrix_cl<double> A_cl(A);
    stan::math::matrix_cl<double> eigenvals_cl
        = stan::math::eigenvalues_sym(A_cl);
    stan::math::matrix_cl<double> eigenvecs_cl
        = stan::math::eigenvectors_sym(A_cl);

    Eigen::VectorXd eigenvals
        = stan::math::from_matrix_cl<Eigen::VectorXd>(eigenvals_cl);
    Eigen::MatrixXd eigenvecs = stan::math::from_matrix_cl(eigenvecs_cl);
    EXPECT_NEAR_REL(A.diagonal().sum(), eigenvals.sum());
    EXPECT_MATRIX_NEAR(eigenvecs * eigenvecs.transpose(),
                       Eigen::MatrixXd::Identity(size, size), 1e-12);
    EXPECT_MATRIX_NEAR(A * eigenvecs, eigenvecs * eigenvals.asDiagonal(),
                       1e-12);
  }
}

TEST(MathMatrix, symmetric_eigensolver_prim_large) {
  int size = 2000;
  Eigen::MatrixXd A = Eigen::MatrixXd::Random(size, size);
  A += A.transpose().eval();

  stan::math::matrix_cl<double> A_cl(A);
  stan::math::matrix_cl<double> eigenvals_cl
      = stan::math::eigenvalues_sym(A_cl);
  stan::math::matrix_cl<double> eigenvecs_cl
      = stan::math::eigenvectors_sym(A_cl);

  Eigen::VectorXd eigenvals
      = stan::math::from_matrix_cl<Eigen::VectorXd>(eigenvals_cl);
  Eigen::MatrixXd eigenvecs = stan::math::from_matrix_cl(eigenvecs_cl);

  EXPECT_NEAR_REL(A.diagonal().sum(), eigenvals.sum());
  EXPECT_MATRIX_NEAR(eigenvecs * eigenvecs.transpose(),
                     Eigen::MatrixXd::Identity(size, size), 1e-12);
  EXPECT_MATRIX_NEAR(A * eigenvecs, eigenvecs * eigenvals.asDiagonal(), 1e-12);
}

TEST(MathMatrix, symmetric_eigensolver_errors) {
  stan::math::matrix_cl<double> empty;
  EXPECT_THROW(stan::math::eigenvalues_sym(empty), std::invalid_argument);
  EXPECT_THROW(stan::math::eigenvectors_sym(empty), std::invalid_argument);

  stan::math::matrix_cl<double> nonsquare(2, 3);
  EXPECT_THROW(stan::math::eigenvalues_sym(nonsquare), std::invalid_argument);
  EXPECT_THROW(stan::math::eigenvectors_sym(nonsquare), std::invalid_argument);

  Eigen::MatrixXd nonsymmetric(2, 2);
  nonsymmetric << 1, 2, 3, 4;
  stan::math::matrix_cl<double> nonsymmetric_cl(nonsymmetric);
  EXPECT_THROW(stan::math::eigenvalues_sym(nonsymmetric_cl), std::domain_error);
  EXPECT_THROW(stan::math::eigenvectors_sym(nonsymmetric_cl),
               std::domain_error);
}

#endif
