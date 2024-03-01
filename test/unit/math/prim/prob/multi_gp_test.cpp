#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(ProbDistributionsMultiGP, MultiGP) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  Matrix<double, Dynamic, 1> mu(5, 1);
  mu.setZero();

  Matrix<double, Dynamic, Dynamic> y(3, 5);
  y << 2.0, -2.0, 11.0, 4.0, -2.0, 11.0, 2.0, -5.0, 11.0, 0.0, -2.0, 11.0, 2.0,
      -2.0, -11.0;

  Matrix<double, Dynamic, Dynamic> Sigma(5, 5);
  Sigma << 9.0, -3.0, 0.0, 0.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0,
      1.0, 0.0, 0.0, 0.0, 1.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0;

  Matrix<double, Dynamic, 1> w(3, 1);
  w << 1.0, 0.5, 1.5;

  double lp_ref(0);
  for (size_t i = 0; i < 3; i++) {
    Matrix<double, Dynamic, 1> cy(y.row(i).transpose());
    Matrix<double, Dynamic, Dynamic> cSigma((1.0 / w[i]) * Sigma);
    lp_ref += stan::math::multi_normal_lpdf(cy, mu, cSigma);
  }

  EXPECT_FLOAT_EQ(lp_ref, stan::math::multi_gp_lpdf(y, Sigma, w));
}

TEST(ProbDistributionsMultiGP, ErrorSigma) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  Matrix<double, Dynamic, Dynamic> y(3, 5);
  y << 2.0, -2.0, 11.0, 4.0, -2.0, 11.0, 2.0, -5.0, 11.0, 0.0, -2.0, 11.0, 2.0,
      -2.0, -11.0;

  Matrix<double, Dynamic, Dynamic> Sigma(5, 5);
  Sigma << 9.0, -3.0, 0.0, 0.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0,
      1.0, 0.0, 0.0, 0.0, 1.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0;

  Matrix<double, Dynamic, 1> w(3, 1);
  w << 1.0, 0.5, 1.5;

  // non-symmetric
  Sigma(0, 1) = -2.5;
  EXPECT_THROW(stan::math::multi_gp_lpdf(y, Sigma, w), std::domain_error);
  Sigma(0, 1) = Sigma(1, 0);

  // non-spd
  Sigma(0, 0) = -3.0;
  EXPECT_THROW(stan::math::multi_gp_lpdf(y, Sigma, w), std::domain_error);
  Sigma(0, 1) = 9.0;
}

TEST(ProbDistributionsMultiGP, ErrorW) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  Matrix<double, Dynamic, Dynamic> y(3, 5);
  y << 2.0, -2.0, 11.0, 4.0, -2.0, 11.0, 2.0, -5.0, 11.0, 0.0, -2.0, 11.0, 2.0,
      -2.0, -11.0;

  Matrix<double, Dynamic, Dynamic> Sigma(5, 5);
  Sigma << 9.0, -3.0, 0.0, 0.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0,
      1.0, 0.0, 0.0, 0.0, 1.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0;

  Matrix<double, Dynamic, 1> w(3, 1);
  w << 1.0, 0.5, 1.5;

  // negative w
  w(0, 0) = -2.5;
  EXPECT_THROW(stan::math::multi_gp_lpdf(y, Sigma, w), std::domain_error);

  // non-finite values
  w(0, 0) = std::numeric_limits<double>::infinity();
  EXPECT_THROW(stan::math::multi_gp_lpdf(y, Sigma, w), std::domain_error);
  w(0, 0) = -std::numeric_limits<double>::infinity();
  EXPECT_THROW(stan::math::multi_gp_lpdf(y, Sigma, w), std::domain_error);
  w(0, 0) = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(stan::math::multi_gp_lpdf(y, Sigma, w), std::domain_error);
}

TEST(ProbDistributionsMultiGP, ErrorY) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  Matrix<double, Dynamic, Dynamic> y(3, 5);
  y << 2.0, -2.0, 11.0, 4.0, -2.0, 11.0, 2.0, -5.0, 11.0, 0.0, -2.0, 11.0, 2.0,
      -2.0, -11.0;

  Matrix<double, Dynamic, Dynamic> Sigma(5, 5);
  Sigma << 9.0, -3.0, 0.0, 0.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0,
      1.0, 0.0, 0.0, 0.0, 1.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0;

  Matrix<double, Dynamic, 1> w(3, 1);
  w << 1.0, 0.5, 1.5;

  // non-finite values
  y(0, 0) = std::numeric_limits<double>::infinity();
  EXPECT_THROW(stan::math::multi_gp_lpdf(y, Sigma, w), std::domain_error);
  y(0, 0) = -std::numeric_limits<double>::infinity();
  EXPECT_THROW(stan::math::multi_gp_lpdf(y, Sigma, w), std::domain_error);
  y(0, 0) = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(stan::math::multi_gp_lpdf(y, Sigma, w), std::domain_error);
}
