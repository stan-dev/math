#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(ErrorHandlingMatrix, isSpsdMatrixPosDef) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y;

  y.resize(3, 3);
  y << 2, -1, 0, -1, 2, -1, 0, -1, 2;
  EXPECT_TRUE(stan::math::is_spsd_matrix(y));

  y << 1, 2, 3, 2, 1, 2, 3, 2, 1;
  EXPECT_FALSE(stan::math::is_spsd_matrix(y));

  y.setZero();
  EXPECT_TRUE(stan::math::is_spsd_matrix(y));
}

TEST(ErrorHandlingMatrix, isSpsdMatrixZero) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y
      = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Zero(3, 3);
  EXPECT_TRUE(stan::math::is_spsd_matrix(y));
}

TEST(ErrorHandlingMatrix, isSpsdNotSquare) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y
      = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Zero(3, 2);
  EXPECT_FALSE(stan::math::is_spsd_matrix(y));
}

TEST(ErrorHandlingMatrix, isSpsdMatrixPosDef_nan) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y;
  double nan = std::numeric_limits<double>::quiet_NaN();

  y.resize(3, 3);
  y << 2, -1, 0, -1, 2, -1, 0, -1, 2;
  EXPECT_TRUE(stan::math::is_spsd_matrix(y));

  y.setZero();
  EXPECT_TRUE(stan::math::is_spsd_matrix(y));

  for (int i = 0; i < y.size(); i++) {
    y << 2, -1, 0, -1, 2, -1, 0, -1, 2;
    y(i) = nan;
    EXPECT_FALSE(stan::math::is_spsd_matrix(y));

    y.setZero();
    y(i) = nan;
    EXPECT_FALSE(stan::math::is_spsd_matrix(y));
  }
}
