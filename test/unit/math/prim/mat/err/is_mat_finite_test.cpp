#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <string>

using stan::math::is_mat_finite;

// ---------- is_mat_finite: matrix tests ----------
TEST(ErrorHandlingScalar, isMatFinite_Matrix) {
  Eigen::Matrix<double, Eigen::Dynamic, 1> x;

  x.resize(3);
  x << -1, 0, 1;
  EXPECT_TRUE(is_mat_finite(x));

  x.resize(3);
  x << -1, 0, std::numeric_limits<double>::infinity();
  EXPECT_FALSE(is_mat_finite(x));

  x.resize(3);
  x << -1, 0, -std::numeric_limits<double>::infinity();
  EXPECT_FALSE(is_mat_finite(x));

  x.resize(3);
  x << -1, 0, std::numeric_limits<double>::quiet_NaN();
  EXPECT_FALSE(is_mat_finite(x));
}

TEST(ErrorHandlingScalar, isMatFinite_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  Eigen::Matrix<double, Eigen::Dynamic, 1> x_mat(3);
  x_mat << nan, 0, 1;
  EXPECT_FALSE(is_mat_finite(x_mat));

  x_mat << 1, nan, 1;
  EXPECT_FALSE(is_mat_finite(x_mat));

  x_mat << 1, 0, nan;
  EXPECT_FALSE(is_mat_finite(x_mat));
}
