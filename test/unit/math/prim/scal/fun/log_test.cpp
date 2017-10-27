#include <gtest/gtest.h>
#include <stan/math.hpp>

TEST(MathFunctions, logInt) {
  using std::log;
  using stan::math::log;
  EXPECT_FLOAT_EQ(std::log(3), log(3));
  EXPECT_FLOAT_EQ(std::log(3.0), log(3.0));
}
