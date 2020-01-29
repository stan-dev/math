#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, inv_logit) {
  using stan::math::inv_logit;
  EXPECT_FLOAT_EQ(0.5, inv_logit(0.0));
  EXPECT_FLOAT_EQ(1.0 / (1.0 + exp(-5.0)), inv_logit(5.0));
}

TEST(MathFunctions, inv_logit_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_TRUE(std::isnan(stan::math::inv_logit(nan)));
}
