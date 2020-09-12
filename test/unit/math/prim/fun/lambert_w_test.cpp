#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <stdexcept>

TEST(MathFunctions, lambert_w) {
  using stan::math::exp;
  using stan::math::lambert_w0;
  using stan::math::lambert_wm1;

  EXPECT_FLOAT_EQ(-1.0, lambert_w0(-1 / exp(1)));
  EXPECT_FLOAT_EQ(1.7455280027406994, lambert_w0(10.));
  EXPECT_FLOAT_EQ(1.7455280027406994, lambert_w0(10));
  EXPECT_FLOAT_EQ(-1.0, lambert_wm1(-1 / exp(1)));
  EXPECT_FLOAT_EQ(lambert_wm1(-std::numeric_limits<double>::min()),
                  -714.96865723796634);
}

TEST(MathFunctions, lambert_wn1_at_0) {
  EXPECT_TRUE(std::isinf(stan::math::lambert_wm1(0)));
}
