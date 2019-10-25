#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(ErrorHandlingScalar, isPositiveFinite) {
  double x = 1;
  EXPECT_TRUE(stan::math::is_positive_finite(x));

  x = -1;
  EXPECT_FALSE(stan::math::is_positive_finite(x));

  x = 0;
  EXPECT_FALSE(stan::math::is_positive_finite(x));

  x = std::numeric_limits<double>::infinity();
  EXPECT_FALSE(stan::math::is_positive_finite(x));

  x = -std::numeric_limits<double>::infinity();
  EXPECT_FALSE(stan::math::is_positive_finite(x));

  x = std::numeric_limits<double>::quiet_NaN();
  EXPECT_FALSE(stan::math::is_positive_finite(x));
}

TEST(ErrorHandlingScalar, isPositiveFinite_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();
  EXPECT_FALSE(stan::math::is_positive_finite(nan));
}
