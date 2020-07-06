#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <string>

TEST(ErrorHandlingMatrix, isCorrMatrix) {
  using stan::math::is_corr_matrix;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y;
  y.resize(2, 2);

  y << 1, 0, 0, 1;
  EXPECT_TRUE(is_corr_matrix(y));

  y << 10, 0, 0, 10;
  EXPECT_FALSE(is_corr_matrix(y));
}

TEST(ErrorHandlingMatrix, isCorrMatrix_nan) {
  using stan::math::is_corr_matrix;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y;
  y.resize(2, 2);
  double nan = std::numeric_limits<double>::quiet_NaN();

  for (int i = 0; i < y.size(); i++) {
    y << 1, 0, 0, 1;
    y(i) = nan;
    EXPECT_FALSE(is_corr_matrix(y));

    y << 10, 0, 0, 10;
    y(i) = nan;
    EXPECT_FALSE(is_corr_matrix(y));
  }
}
