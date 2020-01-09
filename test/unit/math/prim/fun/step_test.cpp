#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(MathFunctions, step_double) {
  using stan::math::step;
  EXPECT_EQ(1.0, step(3.7));
  EXPECT_EQ(1.0, step(0.0));
  EXPECT_EQ(0.0, step(-2.93));
}

TEST(MathFunctions, step_int) {
  using stan::math::step;
  EXPECT_EQ(1.0, step(static_cast<int>(4)));
  EXPECT_EQ(1.0, step(static_cast<int>(0)));
  EXPECT_EQ(0.0, step(static_cast<int>(-3)));
}

TEST(MathFunctions, step_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();
  EXPECT_EQ(1.0, stan::math::step(nan));
}
