#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(MathFunctions, is_uninitialized) {
  using stan::math::is_uninitialized;
  EXPECT_FALSE(is_uninitialized(1.0));
  double x = 1.0;
  EXPECT_FALSE(is_uninitialized(x));

  EXPECT_FALSE(is_uninitialized(3));
  int y = 3;
  EXPECT_FALSE(is_uninitialized(y));
}

TEST(MathFunctions, is_uninitialized_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_FALSE(stan::math::is_uninitialized(nan));
}
