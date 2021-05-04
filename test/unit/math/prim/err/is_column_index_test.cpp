#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(ErrorHandlingMatrix, isColumnIndexMatrix) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y;
  size_t i;

  i = 2;
  y.resize(3, 3);
  EXPECT_TRUE(stan::math::is_column_index(y, i));

  i = 3;
  EXPECT_TRUE(stan::math::is_column_index(y, i));

  y.resize(3, 2);
  EXPECT_FALSE(stan::math::is_column_index(y, i));

  i = 0;
  EXPECT_FALSE(stan::math::is_column_index(y, i));
}

TEST(ErrorHandlingMatrix, isColumnIndexMatrix_nan) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y;
  double nan = std::numeric_limits<double>::quiet_NaN();
  size_t i;

  i = 2;
  y.resize(3, 3);
  y << nan, nan, nan, nan, nan, nan, nan, nan, nan;
  EXPECT_TRUE(stan::math::is_column_index(y, i));

  i = 3;
  EXPECT_TRUE(stan::math::is_column_index(y, i));

  y.resize(3, 2);
  y << nan, nan, nan, nan, nan, nan;
  EXPECT_FALSE(stan::math::is_column_index(y, i));

  i = 0;
  EXPECT_FALSE(stan::math::is_column_index(y, i));
}
