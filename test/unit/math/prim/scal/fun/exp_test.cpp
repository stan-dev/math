#include <stan/math.hpp>
#include <gtest/gtest.h>

TEST(MathFunctions, expInt) {
  using stan::math::exp;
  using std::exp;
  EXPECT_FLOAT_EQ(std::exp(3), exp(3));
  EXPECT_FLOAT_EQ(std::exp(3.0), exp(3.0));
}
