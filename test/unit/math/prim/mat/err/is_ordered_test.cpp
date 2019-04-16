#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <string>

using stan::math::is_ordered;

TEST(ErrorHandlingMatrix, isOrdered) {
  Eigen::Matrix<double, Eigen::Dynamic, 1> y;
  y.resize(3);

  y << 0, 1, 2;
  EXPECT_TRUE(is_ordered(y));

  y << 0, 10, std::numeric_limits<double>::infinity();
  EXPECT_TRUE(is_ordered(y));

  y << -10, 10, std::numeric_limits<double>::infinity();
  EXPECT_TRUE(is_ordered(y));

  y << -std::numeric_limits<double>::infinity(), 10,
      std::numeric_limits<double>::infinity();
  EXPECT_TRUE(is_ordered(y));

  y << 0, 0, 0;
  EXPECT_FALSE(is_ordered(y));

  y << 0, std::numeric_limits<double>::infinity(),
      std::numeric_limits<double>::infinity();
  EXPECT_FALSE(is_ordered(y));

  y << -1, 3, 2;
  EXPECT_FALSE(is_ordered(y));
}

TEST(ErrorHandlingMatrix, isOrdered_one_indexed_message) {
  std::string message;
  Eigen::Matrix<double, Eigen::Dynamic, 1> y;
  y.resize(3);

  y << 0, 5, 1;
  EXPECT_FALSE(is_ordered(y));
}

TEST(ErrorHandlingMatrix, isOrdered_nan) {
  Eigen::Matrix<double, Eigen::Dynamic, 1> y;
  double nan = std::numeric_limits<double>::quiet_NaN();
  y.resize(3);

  y << 0, 1, 2;
  for (int i = 0; i < y.size(); i++) {
    y[i] = nan;
    EXPECT_FALSE(is_ordered(y));
    y[i] = i;
  }
  for (int i = 0; i < y.size(); i++) {
    y << 0, 10, std::numeric_limits<double>::infinity();
    y[i] = nan;
    EXPECT_FALSE(is_ordered(y));
  }
}
