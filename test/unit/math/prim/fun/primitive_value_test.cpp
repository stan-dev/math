#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, primitiveValue) {
  using stan::math::primitive_value;
  // int
  EXPECT_EQ(5, primitive_value(5));
  // uint
  EXPECT_EQ(5U, primitive_value(5U));
  // long >> int
  EXPECT_EQ(10000000000L, primitive_value(10000000000L));
  // char
  EXPECT_EQ('a', primitive_value('a'));

  // double
  EXPECT_EQ(7.3, primitive_value(7.3));
  // float
  EXPECT_EQ(7.3f, primitive_value(7.3f));
}

TEST(MathFunctions, primiviteValueNaN) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_TRUE(std::isnan(stan::math::primitive_value(nan)));
}
