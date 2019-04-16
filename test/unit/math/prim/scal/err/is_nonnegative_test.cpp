#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(ErrorHandlingScalar, isNonnegative) {
  double x = 0;
  EXPECT_TRUE(stan::math::is_nonnegative(x));

  x = std::numeric_limits<double>::infinity();
  EXPECT_TRUE(stan::math::is_nonnegative(x));

  x = -0.01;
  EXPECT_FALSE(stan::math::is_nonnegative(x));

  x = -std::numeric_limits<double>::infinity();
  EXPECT_FALSE(stan::math::is_nonnegative(x));

  x = std::numeric_limits<double>::quiet_NaN();
  EXPECT_FALSE(stan::math::is_nonnegative(x));
}

TEST(ErrorHandlingScalar, isNonnegative_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();
  EXPECT_FALSE(stan::math::is_nonnegative(nan));
}
