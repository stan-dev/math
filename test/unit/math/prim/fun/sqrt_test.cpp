#include <stan/math.hpp>
#include <gtest/gtest.h>

TEST(MathFunctions, sqrtInt) {
  using stan::math::sqrt;
  EXPECT_FLOAT_EQ(std::sqrt(3.0), sqrt(3));
  EXPECT_FLOAT_EQ(std::sqrt(3.1), sqrt(3.1));
  EXPECT_TRUE(stan::math::is_nan(sqrt(-2)));
}
