#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <limits>

TEST(ErrorHandlingMatrix, isSimplex) {
  Eigen::Matrix<double, Eigen::Dynamic, 1> y(2);
  y.setZero();
  y << 0.5, 0.5;

  EXPECT_TRUE(stan::math::is_simplex(y));

  y[1] = 0.55;
  EXPECT_FALSE(stan::math::is_simplex(y));
}

TEST(ErrorHandlingMatrix, isSimplex_negative_value) {
  Eigen::Matrix<double, Eigen::Dynamic, 1> y(100);
  y.setZero();

  y[0] = -0.1;
  y[1] = 1.1;
  EXPECT_FALSE(stan::math::is_simplex(y));

  y.setZero();
  y[0] = 0.1;
  y[1] = -0.1;
  y[2] = 1.0;
  EXPECT_FALSE(stan::math::is_simplex(y));
}

TEST(ErrorHandlingMatrix, isSimplex_sum) {
  Eigen::Matrix<double, Eigen::Dynamic, 1> y(100);
  y.setZero();

  y[13] = 0.9;
  EXPECT_FALSE(stan::math::is_simplex(y));
}

TEST(ErrorHandlingMatrix, isSimplex_length) {
  Eigen::Matrix<double, Eigen::Dynamic, 1> y;

  y.resize(0);
  EXPECT_FALSE(stan::math::is_simplex(y));
}

TEST(ErrorHandlingMatrix, isSimplex_nan) {
  Eigen::Matrix<double, Eigen::Dynamic, 1> y(2);
  y.setZero();
  double nan = std::numeric_limits<double>::quiet_NaN();

  y << nan, 0.5;
  EXPECT_FALSE(stan::math::is_simplex(y));

  y[1] = 0.55;
  EXPECT_FALSE(stan::math::is_simplex(y));

  y << 0.5, nan;
  EXPECT_FALSE(stan::math::is_simplex(y));

  y[0] = nan;
  EXPECT_FALSE(stan::math::is_simplex(y));
}
