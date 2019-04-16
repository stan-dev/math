#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(ErrorHandlingMatrix, isRowIndexMatrix) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y;
  size_t i;

  i = 2;
  y.resize(3, 3);
  EXPECT_TRUE(stan::math::is_row_index(y, i));

  i = 3;
  EXPECT_TRUE(stan::math::is_row_index(y, i));

  y.resize(2, 3);
  EXPECT_FALSE(stan::math::is_row_index(y, i));

  i = 0;
  EXPECT_FALSE(stan::math::is_row_index(y, i));
}

TEST(ErrorHandlingMatrix, isRowIndexMatrix_nan) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y;
  size_t i;
  double nan = std::numeric_limits<double>::quiet_NaN();

  i = 2;
  y.resize(3, 3);
  y << nan, nan, nan, nan, nan, nan, nan, nan, nan;
  EXPECT_TRUE(stan::math::is_row_index(y, i));

  i = 3;
  EXPECT_TRUE(stan::math::is_row_index(y, i));

  y.resize(2, 3);
  y << nan, nan, nan, nan, nan, nan;
  EXPECT_FALSE(stan::math::is_row_index(y, i));

  i = 0;
  EXPECT_FALSE(stan::math::is_row_index(y, i));
}
