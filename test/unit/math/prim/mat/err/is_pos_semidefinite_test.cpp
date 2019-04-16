#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <limits>

using stan::math::is_pos_semidefinite;

TEST(ErrorHandlingMatrix, isPosSemidefinite_size_1) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y;

  y.resize(1, 1);
  y << 0.0;
  EXPECT_TRUE(is_pos_semidefinite(y));

  y << -1.0;
  EXPECT_FALSE(is_pos_semidefinite(y));
}

TEST(ErrorHandlingMatrix, isPosSemidefinite_bad_sizes) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y;

  y.resize(0, 0);
  EXPECT_FALSE(is_pos_semidefinite(y));

  y.resize(2, 3);
  EXPECT_FALSE(is_pos_semidefinite(y));
}

TEST(ErrorHandlingMatrix, isPosSemidefinite) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y;

  y.resize(2, 2);
  y << 1, 0, 0, 1;
  EXPECT_TRUE(is_pos_semidefinite(y));

  y << -1, 0, 0, 1;
  EXPECT_FALSE(is_pos_semidefinite(y));
}

TEST(ErrorHandlingMatrix, isPosSemidefinite_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y;

  y.resize(3, 3);
  y << 2, -1, 0, -1, 2, -1, 0, -1, 2;
  EXPECT_TRUE(is_pos_semidefinite(y));
  for (int i = 0; i < y.rows(); i++)
    for (int j = 0; j < y.cols(); j++) {
      y << 2, -1, 0, -1, 2, -1, 0, -1, 2;
      y(i, j) = nan;
      if (i >= j) {
        EXPECT_FALSE(is_pos_semidefinite(y));
      }
    }
}
