#include <stan/math/prim/scal.hpp>
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

TEST(MathFunctions, is_nan_variadic) {
  using stan::math::is_nan;
  double infinity = std::numeric_limits<double>::infinity();
  double nan = std::numeric_limits<double>::quiet_NaN();
  double min = std::numeric_limits<double>::min();
  double max = std::numeric_limits<double>::max();

  EXPECT_TRUE(is_nan(nan, infinity, 1.0));
  EXPECT_TRUE(is_nan(max, infinity, nan, 2.0, 5.0));
  EXPECT_TRUE(is_nan(max, min, nan));
  EXPECT_FALSE(is_nan(1.0, 2.0, 20.0, 1, 2));
}
