#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(ErrorHandlingScalar, isPositive) {
  using stan::math::is_positive;
  EXPECT_TRUE(is_positive(3.0));
}

TEST(ErrorHandlingScalar, isPositive_nan) {
  using stan::math::is_positive;
  double nan = std::numeric_limits<double>::quiet_NaN();
  EXPECT_FALSE(is_positive(nan));
}
