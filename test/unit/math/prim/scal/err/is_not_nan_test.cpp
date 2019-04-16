#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <limits>

using stan::math::is_not_nan;

TEST(ErrorHandlingScalar, isNotNan) {
  double x = 0;

  EXPECT_TRUE(is_not_nan(x));

  x = std::numeric_limits<double>::infinity();
  EXPECT_TRUE(is_not_nan(x));

  x = -std::numeric_limits<double>::infinity();
  EXPECT_TRUE(is_not_nan(x));

  x = std::numeric_limits<double>::quiet_NaN();
  EXPECT_FALSE(is_not_nan(x));
}
