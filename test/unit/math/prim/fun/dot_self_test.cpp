#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <vector>

TEST(MathFunctions, dot_self) {
  std::vector<double> x(3), y(3);
  x[0] = 2.33;
  x[1] = 8.88;
  x[2] = 9.81;
  y[0] = 2.46;
  y[1] = 4.45;
  y[2] = 1.03;

  EXPECT_FLOAT_EQ(180.5194, stan::math::dot_self(x));
  EXPECT_FLOAT_EQ(26.915, stan::math::dot_self(y));
}

TEST(MathFunctions, dot_self_nan) {
  std::vector<double> x(3);
  x[0] = 2.33;
  x[1] = 8.88;
  x[2] = 9.81;

  double nan = std::numeric_limits<double>::quiet_NaN();
  x[2] = nan;

  EXPECT_TRUE(std::isnan(stan::math::dot_self(x)));

  x[0] = nan;
  x[1] = nan;
  x[2] = nan;
  EXPECT_TRUE(std::isnan(stan::math::dot_self(x)));
}

TEST(MathMatrixPrimMat, dot_self) {
  using stan::math::dot_self;

  Eigen::Matrix<double, Eigen::Dynamic, 1> v1(1);
  v1 << 2.0;
  EXPECT_NEAR(4.0, dot_self(v1), 1E-12);
  Eigen::Matrix<double, Eigen::Dynamic, 1> v2(2);
  v2 << 2.0, 3.0;
  EXPECT_NEAR(13.0, dot_self(v2), 1E-12);
  Eigen::Matrix<double, Eigen::Dynamic, 1> v3(3);
  v3 << 2.0, 3.0, 4.0;
  EXPECT_NEAR(29.0, dot_self(v3), 1E-12);

  Eigen::Matrix<double, 1, Eigen::Dynamic> rv1(1);
  rv1 << 2.0;
  EXPECT_NEAR(4.0, dot_self(rv1), 1E-12);
  Eigen::Matrix<double, 1, Eigen::Dynamic> rv2(2);
  rv2 << 2.0, 3.0;
  EXPECT_NEAR(13.0, dot_self(rv2), 1E-12);
  Eigen::Matrix<double, 1, Eigen::Dynamic> rv3(3);
  rv3 << 2.0, 3.0, 4.0;
  EXPECT_NEAR(29.0, dot_self(rv3), 1E-12);

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> m1(1, 1);
  m1 << 2.0;
  EXPECT_NEAR(4.0, dot_self(m1), 1E-12);
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> m2(2, 1);
  m2 << 2.0, 3.0;
  EXPECT_NEAR(13.0, dot_self(m2), 1E-12);
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> m3(3, 1);
  m3 << 2.0, 3.0, 4.0;
  EXPECT_NEAR(29.0, dot_self(m3), 1E-12);

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> mm2(1, 2);
  mm2 << 2.0, 3.0;
  EXPECT_NEAR(13.0, dot_self(mm2), 1E-12);
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> mm3(1, 3);
  mm3 << 2.0, 3.0, 4.0;
  EXPECT_NEAR(29.0, dot_self(mm3), 1E-12);
}
