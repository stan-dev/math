#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, trunc) {
  using stan::math::trunc;
  EXPECT_FLOAT_EQ(-27, trunc(-27.8239));
  EXPECT_FLOAT_EQ(-1, trunc(-1.5));
  EXPECT_FLOAT_EQ(0, trunc(0));
  EXPECT_FLOAT_EQ(0, trunc(0.0));
  EXPECT_FLOAT_EQ(0, trunc(0.5));
  EXPECT_FLOAT_EQ(27, trunc(27.3239));
}

TEST(MathFunctions, truncNaN) {
  double nan = std::numeric_limits<double>::quiet_NaN();
  EXPECT_TRUE(std::isnan(stan::math::trunc(nan)));
}
