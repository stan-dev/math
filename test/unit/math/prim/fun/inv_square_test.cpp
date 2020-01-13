#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, inv_square) {
  double y = 2.0;
  EXPECT_FLOAT_EQ(1 / (y * y), stan::math::inv_square(y));

  y = 0.0;
  EXPECT_FLOAT_EQ(stan::math::positive_infinity(), stan::math::inv_square(y));

  y = -32.7;
  EXPECT_FLOAT_EQ(1 / (y * y), stan::math::inv_square(y));
}

TEST(MathFunctions, inv_square_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_TRUE(std::isnan(stan::math::inv_square(nan)));
}
