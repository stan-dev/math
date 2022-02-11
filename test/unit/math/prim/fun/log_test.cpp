#include <stan/math.hpp>
#include <gtest/gtest.h>

TEST(MathFunctions, logInt) {
  using stan::math::log;
  EXPECT_FLOAT_EQ(std::log(3), log(3));
  EXPECT_FLOAT_EQ(std::log(3.1), log(3.1));
  EXPECT_FLOAT_EQ(std::log(3.0), log(3.0));
}
