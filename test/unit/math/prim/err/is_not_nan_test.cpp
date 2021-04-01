#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(ErrorHandlingScalar, isNotNan) {
  using stan::math::is_not_nan;
  double x = 0;

  EXPECT_TRUE(is_not_nan(x));

  x = std::numeric_limits<double>::infinity();
  EXPECT_TRUE(is_not_nan(x));

  x = -std::numeric_limits<double>::infinity();
  EXPECT_TRUE(is_not_nan(x));

  x = std::numeric_limits<double>::quiet_NaN();
  EXPECT_FALSE(is_not_nan(x));
}
