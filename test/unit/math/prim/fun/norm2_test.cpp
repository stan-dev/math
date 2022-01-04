#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <vector>

TEST(MathFunctions, norm2) {
  std::vector<double> x(2), y(4);
  x[0] = 3;
  x[1] = 4;

  y[0] = 1;
  y[1] = 2;
  y[2] = 4;
  y[3] = 2;

  EXPECT_FLOAT_EQ(5.0, stan::math::norm2(x));
  EXPECT_FLOAT_EQ(5.0, stan::math::norm2(y));
}

TEST(MathFunctions, norm2_nan) {
  std::vector<double> x(3);
  x[0] = 2.33;
  x[1] = 8.88;
  x[2] = 9.81;

  double nan = std::numeric_limits<double>::quiet_NaN();
  x[2] = nan;

  EXPECT_TRUE(std::isnan(stan::math::norm2(x)));

  x[0] = nan;
  x[1] = nan;
  x[2] = nan;
  EXPECT_TRUE(std::isnan(stan::math::norm2(x)));
}

TEST(MathMatrixPrimMat, norm2) {
  using stan::math::norm2;

  Eigen::Matrix<double, Eigen::Dynamic, 1> v1(1);
  v1 << 2.0;
  EXPECT_NEAR(2.0, norm2(v1), 1E-12);
  Eigen::Matrix<double, Eigen::Dynamic, 1> v2(2);
  v2 << 3.0, 4.0;
  EXPECT_NEAR(5.0, norm2(v2), 1E-12);
  Eigen::Matrix<double, Eigen::Dynamic, 1> v3(4);
  v3 << 1.0, 2.0, 2.0, 4.0;
  EXPECT_NEAR(5.0, norm2(v3), 1E-12);

  Eigen::Matrix<double, 1, Eigen::Dynamic> rv1(1);
  rv1 << 2.0;
  EXPECT_NEAR(2.0, norm2(rv1), 1E-12);
  Eigen::Matrix<double, 1, Eigen::Dynamic> rv2(2);
  rv2 << 3.0, 4.0;
  EXPECT_NEAR(5.0, norm2(rv2), 1E-12);
  Eigen::Matrix<double, 1, Eigen::Dynamic> rv3(4);
  rv3 << 1.0, 2.0, 2.0, 4.0;
  EXPECT_NEAR(5.0, norm2(rv2), 1E-12);

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> m1(1, 1);
  m1 << 2.0;
  EXPECT_NEAR(2.0, norm2(m1), 1E-12);
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> m2(2, 1);
  m2 << 3.0, 4.0;
  EXPECT_NEAR(5.0, norm2(m2), 1E-12);
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> m3(4, 1);
  m3 << 1.0, 2.0, 2.0, 4.0;
  EXPECT_NEAR(5.0, norm2(m3), 1E-12);

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> mm2(1, 2);
  mm2 << 3.0, 4.0;
  EXPECT_NEAR(5.0, norm2(mm2), 1E-12);
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> mm3(1, 4);
  mm3 << 1.0, 2.0, 2.0, 4.0;
  EXPECT_NEAR(5.0, norm2(mm3), 1E-12);
}
