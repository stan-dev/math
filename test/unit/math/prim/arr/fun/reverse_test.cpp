#include <stan/math/prim/arr.hpp>
#include <gtest/gtest.h>

TEST(MathFunctions, reverse) {
  std::vector<double> x(3);
  x[0] = 9;
  x[1] = 2.2;
  x[2] = 3.3;
  std::vector<double> y = stan::math::reverse(x);
  EXPECT_FLOAT_EQ(y[0], 3.3);
  EXPECT_FLOAT_EQ(y[1], 2.2);
  EXPECT_FLOAT_EQ(y[2], 9);
}
