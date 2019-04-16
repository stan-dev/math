#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(ErrorHandlingScalar, isBounded_x) {
  double x = 0;
  double low = -1;
  double high = 1;
  EXPECT_TRUE(stan::math::is_bounded(x, low, high));

  x = low;
  EXPECT_TRUE(stan::math::is_bounded(x, low, high));

  x = high;
  EXPECT_TRUE(stan::math::is_bounded(x, low, high));

  x = low - 1;
  EXPECT_FALSE(stan::math::is_bounded(x, low, high));

  x = high + 1;
  EXPECT_FALSE(stan::math::is_bounded(x, low, high));

  x = std::numeric_limits<double>::quiet_NaN();
  EXPECT_FALSE(stan::math::is_bounded(x, low, high));

  x = -std::numeric_limits<double>::infinity();
  EXPECT_FALSE(stan::math::is_bounded(x, low, high));

  x = std::numeric_limits<double>::infinity();
  EXPECT_FALSE(stan::math::is_bounded(x, low, high));
}

TEST(ErrorHandlingScalar, isBounded_Low) {
  double x = 0;
  double low = -1;
  double high = 1;
  EXPECT_TRUE(stan::math::is_bounded(x, low, high));

  low = -std::numeric_limits<double>::infinity();
  EXPECT_TRUE(stan::math::is_bounded(x, low, high));

  low = std::numeric_limits<double>::quiet_NaN();
  EXPECT_FALSE(stan::math::is_bounded(x, low, high));

  low = std::numeric_limits<double>::infinity();
  EXPECT_FALSE(stan::math::is_bounded(x, low, high));
}

TEST(ErrorHandlingScalar, isBounded_High) {
  double x = 0;
  double low = -1;
  double high = 1;
  EXPECT_TRUE(stan::math::is_bounded(x, low, high));

  high = std::numeric_limits<double>::infinity();
  EXPECT_TRUE(stan::math::is_bounded(x, low, high));

  high = std::numeric_limits<double>::quiet_NaN();
  EXPECT_FALSE(stan::math::is_bounded(x, low, high));

  high = -std::numeric_limits<double>::infinity();
  EXPECT_FALSE(stan::math::is_bounded(x, low, high));
}

TEST(ErrorHandlingScalar, isBounded_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  double x = 0;
  double low = -1;
  double high = 1;
  EXPECT_FALSE(stan::math::is_bounded(nan, low, high));
  EXPECT_FALSE(stan::math::is_bounded(x, nan, high));
  EXPECT_FALSE(stan::math::is_bounded(x, low, nan));
  EXPECT_FALSE(stan::math::is_bounded(nan, nan, high));
  EXPECT_FALSE(stan::math::is_bounded(nan, low, nan));
  EXPECT_FALSE(stan::math::is_bounded(x, nan, nan));
  EXPECT_FALSE(stan::math::is_bounded(nan, nan, nan));
}
