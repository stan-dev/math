#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, fmaxFinite) {
  using stan::math::fmax;
  EXPECT_FLOAT_EQ(1.0, fmax(1, 0));
  EXPECT_FLOAT_EQ(1.0, fmax(1.0, 0));
  EXPECT_FLOAT_EQ(1.0, fmax(1, 0.0));
  EXPECT_FLOAT_EQ(1.0, fmax(1.0, 0.0));

  EXPECT_FLOAT_EQ(1.0, fmax(0, 1));
  EXPECT_FLOAT_EQ(1.0, fmax(0, 1.0));
  EXPECT_FLOAT_EQ(1.0, fmax(0.0, 1));
  EXPECT_FLOAT_EQ(1.0, fmax(0.0, 1.0));
}

TEST(MathFunctions, fmaxNaN) {
  using stan::math::fmax;
  double nan = std::numeric_limits<double>::quiet_NaN();
  EXPECT_FLOAT_EQ(1.0, fmax(1, nan));
  EXPECT_FLOAT_EQ(1.0, fmax(nan, 1));
  EXPECT_TRUE(std::isnan(stan::math::fmax(nan, nan)));
}

TEST(MathFunctions, fmaxInf) {
  using stan::math::fmax;
  double inf = std::numeric_limits<double>::infinity();
  EXPECT_FLOAT_EQ(inf, fmax(inf, 1));
  EXPECT_FLOAT_EQ(inf, fmax(1, inf));
  EXPECT_FLOAT_EQ(inf, fmax(inf, -inf));
  EXPECT_FLOAT_EQ(inf, fmax(-inf, inf));
}
