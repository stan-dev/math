#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(MathFunctions, is_nan) {
  using stan::math::is_nan;
  double infinity = std::numeric_limits<double>::infinity();
  double nan = std::numeric_limits<double>::quiet_NaN();
  double min = std::numeric_limits<double>::min();
  double max = std::numeric_limits<double>::max();
  EXPECT_TRUE(is_nan(nan));
  EXPECT_FALSE(is_nan(infinity));
  EXPECT_FALSE(is_nan(0));
  EXPECT_FALSE(is_nan(1));
  EXPECT_FALSE(is_nan(min));
  EXPECT_FALSE(is_nan(max));
}
