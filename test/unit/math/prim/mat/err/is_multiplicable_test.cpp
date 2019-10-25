#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(ErrorHandlingMatrix, isMultiplicableMatrix) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x;

  y.resize(3, 3);
  x.resize(3, 3);
  EXPECT_TRUE(stan::math::is_multiplicable(x, y));

  x.resize(3, 2);
  y.resize(2, 4);
  EXPECT_TRUE(stan::math::is_multiplicable(x, y));

  y.resize(1, 2);
  EXPECT_FALSE(stan::math::is_multiplicable(x, y));

  x.resize(2, 2);
  EXPECT_FALSE(stan::math::is_multiplicable(x, y));
}

TEST(ErrorHandlingMatrix, isMultiplicableMatrix_0) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x;

  x.resize(3, 0);
  y.resize(0, 3);
  EXPECT_FALSE(stan::math::is_multiplicable(x, y));

  x.resize(0, 4);
  y.resize(4, 3);
  EXPECT_FALSE(stan::math::is_multiplicable(x, y));

  x.resize(3, 4);
  y.resize(4, 0);
  EXPECT_FALSE(stan::math::is_multiplicable(x, y));
}

TEST(ErrorHandlingMatrix, isMultiplicableMatrix_nan) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x;
  double nan = std::numeric_limits<double>::quiet_NaN();

  y.resize(3, 3);
  x.resize(3, 3);
  y << nan, nan, nan, nan, nan, nan, nan, nan, nan;
  x << nan, nan, nan, nan, nan, nan, nan, nan, nan;
  EXPECT_TRUE(stan::math::is_multiplicable(x, y));

  x.resize(3, 2);
  y.resize(2, 4);
  y << nan, nan, nan, nan, nan, nan, nan, nan;
  x << nan, nan, nan, nan, nan, nan;
  EXPECT_TRUE(stan::math::is_multiplicable(x, y));

  y.resize(1, 2);
  y << nan, nan;
  EXPECT_FALSE(stan::math::is_multiplicable(x, y));

  x.resize(2, 2);
  x << nan, nan, nan, nan;
  EXPECT_FALSE(stan::math::is_multiplicable(x, y));
}
