#include <stan/math/prim/scal.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(MathFunctions, step_double) {
  using stan::math::step;
  EXPECT_EQ(0.0, step(-1.0));
  EXPECT_EQ(0.0, step(0.0));
  EXPECT_EQ(1.0, step(0.00000000001));
  EXPECT_EQ(1.0, step(100.0));
}

TEST(MathFunctions, step_int) {
  using stan::math::step;
  EXPECT_EQ(0.0, step(static_cast<int>(-1)));
  EXPECT_EQ(0.0, step(static_cast<int>(0)));
  EXPECT_EQ(1.0, step(static_cast<int>(100)));
}

TEST(MathFunctions, step_inf) {
  using stan::math::step;
  EXPECT_EQ(1.0, step(std::numeric_limits<double>::infinity()));
  EXPECT_EQ(0.0, step(-std::numeric_limits<double>::infinity()));
}

TEST(MathFunctions, step_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();
  EXPECT_EQ(0.0, stan::math::step(nan));
}
