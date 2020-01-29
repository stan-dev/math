#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, logit) {
  using stan::math::logit;
  EXPECT_FLOAT_EQ(0.0, logit(0.5));
  EXPECT_FLOAT_EQ(5.0, logit(1.0 / (1.0 + exp(-5.0))));
}

TEST(MathFunctions, logit_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_TRUE(std::isnan(stan::math::logit(nan)));
}
