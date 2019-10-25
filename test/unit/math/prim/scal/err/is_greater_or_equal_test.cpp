#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(ErrorHandlingScalar, isGreaterOrEqual) {
  double x = 10.0;
  double lb = 0.0;

  EXPECT_TRUE(stan::math::is_greater_or_equal(x, lb));

  x = -1.0;
  EXPECT_FALSE(stan::math::is_greater_or_equal(x, lb));

  x = lb;
  EXPECT_TRUE(stan::math::is_greater_or_equal(x, lb));

  x = std::numeric_limits<double>::infinity();
  EXPECT_TRUE(stan::math::is_greater_or_equal(x, lb));

  x = 10.0;
  lb = std::numeric_limits<double>::infinity();
  EXPECT_FALSE(stan::math::is_greater_or_equal(x, lb));

  x = std::numeric_limits<double>::infinity();
  lb = std::numeric_limits<double>::infinity();
  EXPECT_TRUE(stan::math::is_greater_or_equal(x, lb));
}

TEST(ErrorHandlingScalar, isGreaterOrEqual_nan) {
  double x = 10.0;
  double lb = 0.0;
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_FALSE(stan::math::is_greater_or_equal(nan, lb));
  EXPECT_FALSE(stan::math::is_greater_or_equal(x, nan));
  EXPECT_FALSE(stan::math::is_greater_or_equal(nan, nan));
}
