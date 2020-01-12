#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, log1p_exp) {
  using stan::math::log1p_exp;

  // exp(10000.0) overflows
  EXPECT_FLOAT_EQ(10000.0, log1p_exp(10000.0));
  EXPECT_FLOAT_EQ(0.0, log1p_exp(-10000.0));
}

TEST(MathFunctions, log1p_exp_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_TRUE(std::isnan(stan::math::log1p_exp(nan)));
}
