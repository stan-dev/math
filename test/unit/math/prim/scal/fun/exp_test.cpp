#include <gtest/gtest.h>
#include <stan/math.hpp>

TEST(MathFunctions, expInt) {
  using stan::math::exp;
  using std::exp;
  EXPECT_FLOAT_EQ(std::exp(3), exp(3));
  EXPECT_FLOAT_EQ(std::exp(3.0), exp(3.0));
}
