#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(ErrorHandlingScalar, isNotNan) {
  double x = 0;

  EXPECT_TRUE(stan::math::is_not_nan(x));

  x = std::numeric_limits<double>::infinity();
  EXPECT_TRUE(stan::math::is_not_nan(x));

  x = -std::numeric_limits<double>::infinity();
  EXPECT_TRUE(stan::math::is_not_nan(x));

  x = std::numeric_limits<double>::quiet_NaN();
  EXPECT_FALSE(stan::math::is_not_nan(x));
}
