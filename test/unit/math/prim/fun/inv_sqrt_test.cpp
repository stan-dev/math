#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, inv_sqrt) {
  double y = 4.0;
  EXPECT_FLOAT_EQ(1 / 2.0, stan::math::inv_sqrt(y));

  y = 25.0;
  EXPECT_FLOAT_EQ(1 / 5.0, stan::math::inv_sqrt(y));

  y = 0.0;
  EXPECT_FLOAT_EQ(stan::math::positive_infinity(), stan::math::inv_sqrt(y));

  y = -50.0;
  EXPECT_TRUE(std::isnan(stan::math::inv_sqrt(y)));
}

TEST(MathFunctions, inv_sqrt_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_TRUE(std::isnan(stan::math::inv_sqrt(nan)));
}
