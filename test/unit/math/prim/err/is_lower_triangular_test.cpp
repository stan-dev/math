#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <string>

TEST(ErrorHandlingMatrix, isLowerTriangular) {
  using stan::math::is_lower_triangular;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y;

  y.resize(1, 1);
  y << 1;
  EXPECT_TRUE(is_lower_triangular(y));

  y.resize(1, 2);
  y << 1, 0;
  EXPECT_TRUE(is_lower_triangular(y));

  y(0, 1) = 1;
  EXPECT_FALSE(is_lower_triangular(y));

  y.resize(2, 2);
  y << 1, 0, 2, 3;
  EXPECT_TRUE(is_lower_triangular(y));

  y << 1, 2, 3, 4;
  EXPECT_FALSE(is_lower_triangular(y));

  y.resize(3, 2);
  y << 1, 0, 2, 3, 4, 5;
  EXPECT_TRUE(is_lower_triangular(y));

  y(0, 1) = 1.5;
  EXPECT_FALSE(is_lower_triangular(y));

  y.resize(2, 3);
  y << 1, 0, 0, 4, 5, 0;
  EXPECT_TRUE(is_lower_triangular(y));
}

TEST(ErrorHandlingMatrix, isLowerTriangular_nan) {
  using stan::math::is_lower_triangular;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y;
  double nan = std::numeric_limits<double>::quiet_NaN();

  y.resize(1, 1);
  y << nan;
  EXPECT_TRUE(is_lower_triangular(y));

  y.resize(1, 2);
  y << nan, 0;
  EXPECT_TRUE(is_lower_triangular(y));

  y(0, 1) = nan;
  EXPECT_FALSE(is_lower_triangular(y));

  y.resize(2, 2);
  y << nan, 0, nan, nan;
  EXPECT_TRUE(is_lower_triangular(y));

  y << 1, nan, nan, 4;
  EXPECT_FALSE(is_lower_triangular(y));

  y.resize(3, 2);
  y << nan, 0, 2, nan, 4, 5;
  EXPECT_TRUE(is_lower_triangular(y));

  y(0, 1) = nan;
  EXPECT_FALSE(is_lower_triangular(y));

  y.resize(2, 3);
  y << nan, 0, 0, 4, nan, 0;
  EXPECT_TRUE(is_lower_triangular(y));

  y(0, 2) = nan;
  EXPECT_FALSE(is_lower_triangular(y));
}
