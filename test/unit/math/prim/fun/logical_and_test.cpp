#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(MathFunctions, logical_and) {
  using stan::math::logical_and;
  EXPECT_TRUE(logical_and(1, 1));
  EXPECT_TRUE(logical_and(5.7, -1.9));

  EXPECT_FALSE(logical_and(0, 0));
  EXPECT_FALSE(logical_and(0, 1));
  EXPECT_FALSE(logical_and(1, 0));
  EXPECT_FALSE(logical_and(0.0, 0.0));
  EXPECT_FALSE(logical_and(0.0, 1.0));
  EXPECT_FALSE(logical_and(1, 0.0));
}

TEST(MathFunctions, logical_and_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_TRUE(stan::math::logical_and(1.0, nan));
  EXPECT_TRUE(stan::math::logical_and(nan, 2.0));
  EXPECT_TRUE(stan::math::logical_and(nan, nan));
}
