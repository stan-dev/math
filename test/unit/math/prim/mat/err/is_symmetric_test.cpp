#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <string>

TEST(ErrorHandlingMatrix, isSymmetric) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y;

  y.resize(2, 2);
  y << 1, 3, 3, 1;
  EXPECT_TRUE(stan::math::is_symmetric(y));

  y << 1, 3.5, 3, 1;
  EXPECT_FALSE(stan::math::is_symmetric(y));
}

TEST(ErrorHandlingMatrix, isSymmetric_one_indexed_message) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y;
  std::string message;

  y.resize(2, 2);
  y << 1, 0, 3, 1;
  EXPECT_FALSE(stan::math::is_symmetric(y));
}

TEST(ErrorHandlingMatrix, isSymmetric_nan) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y;
  double nan = std::numeric_limits<double>::quiet_NaN();

  y.resize(2, 2);
  y << 1, nan, 3, 1;
  EXPECT_FALSE(stan::math::is_symmetric(y));

  y << nan, 3, 3, 1;
  EXPECT_TRUE(stan::math::is_symmetric(y));

  y.resize(1, 1);
  y << nan;
  EXPECT_TRUE(stan::math::is_symmetric(y));
}

TEST(ErrorHandlingMatrix, isSymmetric_non_square) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y;

  y.resize(2, 3);
  EXPECT_FALSE(stan::math::is_symmetric(y));
}
