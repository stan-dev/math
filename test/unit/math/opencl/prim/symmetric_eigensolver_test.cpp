#ifdef STAN_OPENCL
#include <iostream>
#include <stan/math/opencl/mrrr.hpp>
#include <test/unit/math/expect_near_rel.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

//TEST(MathMatrix, symmetric_eigensolver_small) {
//  int size = 7;
//  Eigen::MatrixXd A = Eigen::MatrixXd::Random(size, size);
//  //  Eigen::MatrixXd A = Eigen::MatrixXd::Constant(size, size, 1);
//  A += A.transpose().eval();

//  stan::math::matrix_cl<double> A_cl(A);
//  stan::math::matrix_cl<double> eigenvals_cl;
//  stan::math::matrix_cl<double> eigenvecs_cl;

//  stan::math::symmetric_eigensolver(A_cl, eigenvals_cl, eigenvecs_cl);

//  Eigen::VectorXd eigenvals
//      = stan::math::from_matrix_cl<Eigen::VectorXd>(eigenvals_cl);
//  Eigen::MatrixXd eigenvecs = stan::math::from_matrix_cl(eigenvecs_cl);
//  EXPECT_NEAR_REL(A.diagonal().sum(), eigenvals.sum());
//  EXPECT_MATRIX_NEAR(eigenvecs * eigenvecs.transpose(),
//                     Eigen::MatrixXd::Identity(size, size), 1e-12);
//  EXPECT_MATRIX_NEAR(A * eigenvecs, eigenvecs * eigenvals.asDiagonal(), 1e-12);
//}

TEST(MathMatrix, symmetric_eigensolver_large) {
  int size = 2000;
  Eigen::MatrixXd A = Eigen::MatrixXd::Random(size, size).array();
  A += A.transpose().eval();

  stan::math::matrix_cl<double> A_cl(A);
  A_glob_cl = A_cl;
  stan::math::matrix_cl<double> eigenvals_cl;
  stan::math::matrix_cl<double> eigenvecs_cl;

  stan::math::symmetric_eigensolver(A_cl, eigenvals_cl, eigenvecs_cl);

  std::cout << "\n\nfirst done\n\n";

  Eigen::VectorXd eigenvals
      = stan::math::from_matrix_cl<Eigen::VectorXd>(eigenvals_cl);
  Eigen::MatrixXd eigenvecs = stan::math::from_matrix_cl(eigenvecs_cl);

  stan::math::matrix_cl<double> eigenvals2_cl;
  stan::math::matrix_cl<double> eigenvecs2_cl;

  stan::math::symmetric_eigensolver(A_cl, eigenvals2_cl, eigenvecs2_cl);
  std::cout << "\n\nsecond done\n\n";
//  stan::math::symmetric_eigensolver(A_cl, eigenvals2_cl, eigenvecs2_cl);
//  std::cout << "\n\nthird done\n\n";
//  stan::math::symmetric_eigensolver(A_cl, eigenvals2_cl, eigenvecs2_cl);
//  std::cout << "\n\nfourth done\n\n";
//  stan::math::symmetric_eigensolver(A_cl, eigenvals2_cl, eigenvecs2_cl);
//  std::cout << "\n\nfifth done\n\n";
  Eigen::VectorXd eigenvals2
      = stan::math::from_matrix_cl<Eigen::VectorXd>(eigenvals2_cl);
  Eigen::MatrixXd eigenvecs2 = stan::math::from_matrix_cl(eigenvecs2_cl);
  //  EXPECT_MATRIX_NEAR(stan::math::from_matrix_cl(eigenvals_cl),
  //                     stan::math::from_matrix_cl(eigenvals2_cl), 1e-12);
  //  EXPECT_MATRIX_NEAR(stan::math::from_matrix_cl(eigenvecs_cl),
  //                     stan::math::from_matrix_cl(eigenvecs2_cl), 1e-12);
  std::cout << eigenvals.head(10).transpose() << std::endl;
  std::cout << eigenvals2.head(10).transpose() << std::endl;
  std::cout << eigenvals.tail(10).transpose() << std::endl;
  std::cout << eigenvals2.tail(10).transpose() << std::endl;
  std::cout << "eigenval diff: "
            << (eigenvals - eigenvals2).array().abs().maxCoeff()
            << std::endl;
  std::cout << "eigenval avg diff: "
            << (eigenvals - eigenvals2).array().abs().mean()
            << std::endl;
  std::cout << "eigenvec diff: "
            << (eigenvecs - eigenvecs2).array().abs().maxCoeff()
            << std::endl;
  std::cout << "eigenvec avg diff: "
            << (eigenvecs - eigenvecs2).array().abs().mean()
            << std::endl;

  std::cout << "trace rel err: " << abs(A.diagonal().sum() - eigenvals.sum())
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
      << ((A * eigenvecs - eigenvecs * eigenvals.asDiagonal()).array().abs())
             .array()
             .abs()
             .maxCoeff()
      << std::endl;

  std::cout << "trace2 rel err: " << abs(A.diagonal().sum() - eigenvals2.sum())
            << std::endl;
  std::cout << "vec ortho2 err: "
            << ((eigenvecs2 * eigenvecs2.transpose())
                - Eigen::MatrixXd::Identity(size, size))
                   .array()
                   .abs()
                   .maxCoeff()
            << std::endl;
  std::cout
      << "eigen eq2 err: "
      << ((A * eigenvecs2 - eigenvecs2 * eigenvals2.asDiagonal()).array().abs())
             .array()
             .abs()
             .maxCoeff()
      << std::endl;
//  EXPECT_NEAR_REL(A.diagonal().sum(), eigenvals.sum());
//  EXPECT_MATRIX_NEAR(eigenvecs * eigenvecs.transpose(),
//                     Eigen::MatrixXd::Identity(size, size), 1e-12);
//  EXPECT_MATRIX_NEAR(A * eigenvecs, eigenvecs * eigenvals.asDiagonal(), 1e-12);
}

//TEST(MathMatrix, symmetric_eigensolver_prim_small) {
//  int size = 7;
//  Eigen::MatrixXd A = Eigen::MatrixXd::Random(size, size);
//  A += A.transpose().eval();

//  stan::math::matrix_cl<double> A_cl(A);
//  stan::math::matrix_cl<double> eigenvals_cl;
//  stan::math::matrix_cl<double> eigenvecs_cl;

//  stan::math::symmetric_eigensolver(A_cl, eigenvals_cl, eigenvecs_cl);
//  stan::math::matrix_cl<double> eigenvals2_cl
//      = stan::math::eigenvalues_sym(A_cl);
//  stan::math::matrix_cl<double> eigenvecs2_cl
//      = stan::math::eigenvectors_sym(A_cl);
//  EXPECT_MATRIX_NEAR(stan::math::from_matrix_cl(eigenvals_cl),
//                     stan::math::from_matrix_cl(eigenvals2_cl), 1e-12);
//  EXPECT_MATRIX_NEAR(stan::math::from_matrix_cl(eigenvecs_cl),
//                     stan::math::from_matrix_cl(eigenvecs2_cl), 1e-12);

//  Eigen::VectorXd eigenvals
//      = stan::math::from_matrix_cl<Eigen::VectorXd>(eigenvals_cl);
//  Eigen::MatrixXd eigenvecs = stan::math::from_matrix_cl(eigenvecs_cl);
//  EXPECT_NEAR_REL(A.diagonal().sum(), eigenvals.sum());
//  EXPECT_MATRIX_NEAR(eigenvecs * eigenvecs.transpose(),
//                     Eigen::MatrixXd::Identity(size, size), 1e-12);
//  EXPECT_MATRIX_NEAR(A * eigenvecs, eigenvecs * eigenvals.asDiagonal(), 1e-12);
//}

//TEST(MathMatrix, symmetric_eigensolver_prim_large) {
//  int size = 2000;
//  Eigen::MatrixXd A = Eigen::MatrixXd::Random(size, size);
//  A += A.transpose().eval();

//  stan::math::matrix_cl<double> A_cl(A);
//  stan::math::matrix_cl<double> eigenvals_cl
//      = stan::math::eigenvalues_sym(A_cl);
//  std::cout << "\n\nfirst eigenvals done\n\n";
//  stan::math::matrix_cl<double> eigenvecs_cl
//      = stan::math::eigenvectors_sym(A_cl);
//  std::cout << "\n\nfirst eigenvecs done\n\n";

//  Eigen::VectorXd eigenvals
//      = stan::math::from_matrix_cl<Eigen::VectorXd>(eigenvals_cl);
//  Eigen::MatrixXd eigenvecs = stan::math::from_matrix_cl(eigenvecs_cl);

//  stan::math::matrix_cl<double> eigenvals2_cl
//      = stan::math::eigenvalues_sym(A_cl);
//  std::cout << "\n\nsecond eigenvals done\n\n";
//  stan::math::matrix_cl<double> eigenvecs2_cl
//      = stan::math::eigenvectors_sym(A_cl);
//  std::cout << "\n\nsecond eigenvecs done\n\n";
//  Eigen::VectorXd eigenvals2
//      = stan::math::from_matrix_cl<Eigen::VectorXd>(eigenvals2_cl);
//  Eigen::MatrixXd eigenvecs2 = stan::math::from_matrix_cl(eigenvecs2_cl);
//  //  EXPECT_MATRIX_NEAR(stan::math::from_matrix_cl(eigenvals_cl),
//  //                     stan::math::from_matrix_cl(eigenvals2_cl), 1e-12);
//  //  EXPECT_MATRIX_NEAR(stan::math::from_matrix_cl(eigenvecs_cl),
//  //                     stan::math::from_matrix_cl(eigenvecs2_cl), 1e-12);
//  std::cout << "eigenval diff: "
//            << (eigenvals - eigenvals2).array().abs().maxCoeff()
//            << std::endl;
//  std::cout << "eigenval avg diff: "
//            << (eigenvals - eigenvals2).array().abs().mean()
//            << std::endl;
//  std::cout << "eigenvec diff: "
//            << (eigenvecs - eigenvecs2).array().abs().maxCoeff()
//            << std::endl;
//  std::cout << "eigenvec avg diff: "
//            << (eigenvecs - eigenvecs2).array().abs().mean()
//            << std::endl;

//  std::cout << "trace rel err: " << abs(A.diagonal().sum() - eigenvals.sum())
//            << std::endl;
//  std::cout << "vec ortho err: "
//            << ((eigenvecs * eigenvecs.transpose())
//                - Eigen::MatrixXd::Identity(size, size))
//                   .array()
//                   .abs()
//                   .maxCoeff()
//            << std::endl;
//  std::cout
//      << "eigen eq err: "
//      << ((A * eigenvecs - eigenvecs * eigenvals.asDiagonal()).array().abs())
//             .array()
//             .abs()
//             .maxCoeff()
//      << std::endl;

//  std::cout << "trace2 rel err: " << abs(A.diagonal().sum() - eigenvals2.sum())
//            << std::endl;
//  std::cout << "vec ortho2 err: "
//            << ((eigenvecs2 * eigenvecs2.transpose())
//                - Eigen::MatrixXd::Identity(size, size))
//                   .array()
//                   .abs()
//                   .maxCoeff()
//            << std::endl;
//  std::cout
//      << "eigen eq2 err: "
//      << ((A * eigenvecs2 - eigenvecs2 * eigenvals2.asDiagonal()).array().abs())
//             .array()
//             .abs()
//             .maxCoeff()
//      << std::endl;
////  EXPECT_NEAR_REL(A.diagonal().sum(), eigenvals.sum());
////  EXPECT_MATRIX_NEAR(eigenvecs * eigenvecs.transpose(),
////                     Eigen::MatrixXd::Identity(size, size), 1e-12);
////  EXPECT_MATRIX_NEAR(A * eigenvecs, eigenvecs * eigenvals.asDiagonal(), 1e-12);
//}

#endif
