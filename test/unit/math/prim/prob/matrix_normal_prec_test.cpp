#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(ProbDistributionsMatrixNormal, MatrixNormalPrec) {
  using Eigen::MatrixXd;
  MatrixXd y(3, 5);
  y << 2.0, -2.0, 11.0, 4.0, -2.0, 11.0, 2.0, -5.0, 11.0, 0.0, -2.0, 11.0, 2.0,
      -2.0, -11.0;

  MatrixXd mu = MatrixXd::Zero(3, 5);

  MatrixXd Sigma(3, 3);
  Sigma << 1.0, 0.5, 0.1, 0.5, 1.0, 0.2, 0.1, 0.2, 1.0;

  MatrixXd D(5, 5);
  D << 9.0, -3.0, 0.0, 0.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 1.0,
      0.0, 0.0, 0.0, 1.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0;

  double lp_ref;
  lp_ref = stan::math::matrix_normal_prec_lpdf(y, mu, Sigma, D);
  EXPECT_FLOAT_EQ(lp_ref, -2132.0748232368409845);
}

TEST(ProbDistributionsMatrixNormal, ErrorSigma) {
  using Eigen::MatrixXd;
  MatrixXd y(3, 5);
  y << 2.0, -2.0, 11.0, 4.0, -2.0, 11.0, 2.0, -5.0, 11.0, 0.0, -2.0, 11.0, 2.0,
      -2.0, -11.0;

  MatrixXd mu = MatrixXd::Zero(3, 5);

  MatrixXd Sigma(3, 3);
  Sigma << 1.0, 0.5, 0.1, 0.5, 1.0, 0.2, 0.1, 0.2, 1.0;

  MatrixXd D(5, 5);
  D << 9.0, -3.0, 0.0, 0.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 1.0,
      0.0, 0.0, 0.0, 1.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0;

  // non-symmetric
  Sigma(0, 1) = -2.5;
  EXPECT_THROW(stan::math::matrix_normal_prec_lpdf(y, mu, Sigma, D),
               std::domain_error);
  Sigma(0, 1) = Sigma(1, 0);

  // non-spd
  Sigma(0, 0) = -3.0;
  EXPECT_THROW(stan::math::matrix_normal_prec_lpdf(y, mu, Sigma, D),
               std::domain_error);
  Sigma(0, 0) = 9.0;
}

TEST(ProbDistributionsMatrixNormal, ErrorD) {
  using Eigen::MatrixXd;
  MatrixXd y(3, 5);
  y << 2.0, -2.0, 11.0, 4.0, -2.0, 11.0, 2.0, -5.0, 11.0, 0.0, -2.0, 11.0, 2.0,
      -2.0, -11.0;

  MatrixXd mu = MatrixXd::Zero(3, 5);

  MatrixXd Sigma(3, 3);
  Sigma << 1.0, 0.5, 0.1, 0.5, 1.0, 0.2, 0.1, 0.2, 1.0;

  MatrixXd D(5, 5);
  D << 9.0, -3.0, 0.0, 0.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 1.0,
      0.0, 0.0, 0.0, 1.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0;

  // non-symmetric
  D(0, 1) = -2.5;
  EXPECT_THROW(stan::math::matrix_normal_prec_lpdf(y, mu, Sigma, D),
               std::domain_error);
  D(0, 1) = Sigma(1, 0);

  // non-spd
  D(0, 0) = -3.0;
  EXPECT_THROW(stan::math::matrix_normal_prec_lpdf(y, mu, Sigma, D),
               std::domain_error);
  D(0, 0) = 1.0;
}

TEST(ProbDistributionsMatrixNormal, ErrorY) {
  using Eigen::MatrixXd;
  MatrixXd y(3, 5);
  y << 2.0, -2.0, 11.0, 4.0, -2.0, 11.0, 2.0, -5.0, 11.0, 0.0, -2.0, 11.0, 2.0,
      -2.0, -11.0;

  MatrixXd mu = MatrixXd::Zero(3, 5);

  MatrixXd Sigma(3, 3);
  Sigma << 1.0, 0.5, 0.1, 0.5, 1.0, 0.2, 0.1, 0.2, 1.0;

  MatrixXd D(5, 5);
  D << 9.0, -3.0, 0.0, 0.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 1.0,
      0.0, 0.0, 0.0, 1.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0;

  // non-finite values
  y(0, 0) = std::numeric_limits<double>::infinity();
  EXPECT_THROW(stan::math::matrix_normal_prec_lpdf(y, mu, Sigma, D),
               std::domain_error);
  y(0, 0) = -std::numeric_limits<double>::infinity();
  EXPECT_THROW(stan::math::matrix_normal_prec_lpdf(y, mu, Sigma, D),
               std::domain_error);
  y(0, 0) = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(stan::math::matrix_normal_prec_lpdf(y, mu, Sigma, D),
               std::domain_error);
}
