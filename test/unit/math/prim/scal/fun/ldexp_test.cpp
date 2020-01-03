#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, ldexp_double) {
  using stan::math::ldexp;

  EXPECT_FLOAT_EQ(0.0, ldexp(0.0, 5));
  EXPECT_FLOAT_EQ(32.0, ldexp(1.0, 5));
  EXPECT_FLOAT_EQ(64.0, ldexp(2.0, 5));
  EXPECT_FLOAT_EQ(96.0, ldexp(3.0, 5));

  EXPECT_FLOAT_EQ(-32.0, ldexp(-1.0, 5));
  EXPECT_FLOAT_EQ(-64.0, ldexp(-2.0, 5));
  EXPECT_FLOAT_EQ(-96.0, ldexp(-3.0, 5));
}

TEST(MathFunctions, ldexp_int) {
  using stan::math::ldexp;

  EXPECT_FLOAT_EQ(0.0, ldexp(static_cast<int>(0), 5));
  EXPECT_FLOAT_EQ(32.0, ldexp(static_cast<int>(1), 5));
  EXPECT_FLOAT_EQ(64.0, ldexp(static_cast<int>(2), 5));
  EXPECT_FLOAT_EQ(96.0, ldexp(static_cast<int>(3), 5));

  EXPECT_FLOAT_EQ(-32.0, ldexp(static_cast<int>(-1), 5));
  EXPECT_FLOAT_EQ(-64.0, ldexp(static_cast<int>(-2), 5));
  EXPECT_FLOAT_EQ(-96.0, ldexp(static_cast<int>(-3), 5));
}

TEST(MathFunctions, ldexp_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_TRUE(std::isnan(stan::math::ldexp(nan, 5)));
}
