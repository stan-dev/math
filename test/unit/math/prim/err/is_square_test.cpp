#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(ErrorHandlingMatrix, isSquareMatrix) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y;

  y.resize(3, 3);
  EXPECT_TRUE(stan::math::is_square(y));

  y.resize(3, 2);
  EXPECT_FALSE(stan::math::is_square(y));
}

TEST(ErrorHandlingMatrix, isSquareMatrix_nan) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y;
  double nan = std::numeric_limits<double>::quiet_NaN();

  y.resize(3, 3);
  y << nan, nan, nan, nan, nan, nan, nan, nan, nan;
  EXPECT_TRUE(stan::math::is_square(y));

  y.resize(3, 2);
  y << nan, nan, nan, nan, nan, nan;
  EXPECT_FALSE(stan::math::is_square(y));
}

TEST(ErrorHandlingMatrix, isSquareMatrix_0x0) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y;

  y.resize(0, 0);
  EXPECT_TRUE(stan::math::is_square(y));
}

TEST(ErrorHandlingMatrix, isSquareMatrix_0_size) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y;

  y.resize(0, 10);
  EXPECT_FALSE(stan::math::is_square(y));

  y.resize(10, 0);
  EXPECT_FALSE(stan::math::is_square(y));
}
