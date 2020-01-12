#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, fdim_double) {
  using stan::math::fdim;

  EXPECT_FLOAT_EQ(1.0, fdim(3.0, 2.0));
  EXPECT_FLOAT_EQ(0.0, fdim(2.0, 3.0));

  EXPECT_FLOAT_EQ(2.5, fdim(4.5, 2.0));
}

TEST(MathFunctions, fdim_int) {
  using stan::math::fdim;

  // promotes results to double
  EXPECT_FLOAT_EQ(1.0, fdim(static_cast<int>(3), static_cast<int>(2)));
  EXPECT_FLOAT_EQ(0.0, fdim(static_cast<int>(2), static_cast<int>(3)));
}

TEST(MathFunctions, fdim_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_TRUE(std::isnan(stan::math::fdim(3.0, nan)));

  EXPECT_TRUE(std::isnan(stan::math::fdim(nan, 3.0)));

  EXPECT_TRUE(std::isnan(stan::math::fdim(nan, nan)));
}
