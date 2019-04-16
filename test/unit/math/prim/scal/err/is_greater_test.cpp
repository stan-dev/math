#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(ErrorHandlingScalar, isGreater) {
  double x = 10.0;
  double lb = 0.0;

  EXPECT_TRUE(stan::math::is_greater(x, lb));

  x = -1.0;
  EXPECT_FALSE(stan::math::is_greater(x, lb));

  x = lb;
  EXPECT_FALSE(stan::math::is_greater(x, lb));

  x = std::numeric_limits<double>::infinity();
  EXPECT_TRUE(stan::math::is_greater(x, lb));

  x = 10.0;
  lb = std::numeric_limits<double>::infinity();
  EXPECT_FALSE(stan::math::is_greater(x, lb));

  x = std::numeric_limits<double>::infinity();
  lb = std::numeric_limits<double>::infinity();
  EXPECT_FALSE(stan::math::is_greater(x, lb));
}

TEST(ErrorHandlingScalar, isGreater_nan) {
  double x = 10.0;
  double lb = 0.0;
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_FALSE(stan::math::is_greater(nan, lb));
  EXPECT_FALSE(stan::math::is_greater(x, nan));
  EXPECT_FALSE(stan::math::is_greater(nan, nan));
}
