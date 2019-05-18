#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <limits>

using stan::math::is_less_or_equal;

TEST(ErrorHandlingScalar, isLessOrEqual) {
  double x = -10.0;
  double lb = 0.0;

  EXPECT_TRUE(is_less_or_equal(x, lb));

  x = 1.0;
  EXPECT_FALSE(is_less_or_equal(x, lb));

  x = lb;
  EXPECT_TRUE(is_less_or_equal(x, lb));

  x = -std::numeric_limits<double>::infinity();
  EXPECT_TRUE(is_less_or_equal(x, lb));

  x = -10.0;
  lb = -std::numeric_limits<double>::infinity();
  EXPECT_FALSE(is_less_or_equal(x, lb));

  x = -std::numeric_limits<double>::infinity();
  lb = -std::numeric_limits<double>::infinity();
  EXPECT_TRUE(is_less_or_equal(x, lb));
}

TEST(ErrorHandlingScalar, isLessOrEqual_nan) {
  double x = 10.0;
  double lb = 0.0;
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_FALSE(is_less_or_equal(nan, lb));
  EXPECT_FALSE(is_less_or_equal(x, nan));
  EXPECT_FALSE(is_less_or_equal(nan, nan));
}
