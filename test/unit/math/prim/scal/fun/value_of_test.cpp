#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, value_of) {
  using stan::math::value_of;
  double x = 5.0;
  EXPECT_FLOAT_EQ(5.0, value_of(x));
  EXPECT_FLOAT_EQ(5.0, value_of(5));
}

TEST(MathFunctions, value_of_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_TRUE(std::isnan(stan::math::value_of(nan)));
}
