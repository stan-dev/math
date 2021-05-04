#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(MathFunctions, is_any_nan_variadic_prim) {
  using stan::math::is_any_nan;
  double infinity = std::numeric_limits<double>::infinity();
  double nan = std::numeric_limits<double>::quiet_NaN();
  double min = std::numeric_limits<double>::min();
  double max = std::numeric_limits<double>::max();

  EXPECT_TRUE(is_any_nan(nan, infinity, 1.0));
  EXPECT_TRUE(is_any_nan(max, infinity, nan, 2.0, 5.0));
  EXPECT_TRUE(is_any_nan(max, min, nan));
  EXPECT_FALSE(is_any_nan(1.0, 2.0, 20.0, 1, 2));
}
