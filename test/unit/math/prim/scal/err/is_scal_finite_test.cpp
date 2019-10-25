#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(ErrorHandlingScalar, isScalFinite) {
  double x = 0;
  EXPECT_TRUE(stan::math::is_scal_finite(x));

  x = std::numeric_limits<double>::infinity();
  EXPECT_FALSE(stan::math::is_scal_finite(x));

  x = -std::numeric_limits<double>::infinity();
  EXPECT_FALSE(stan::math::is_scal_finite(x));

  x = std::numeric_limits<double>::quiet_NaN();
  EXPECT_FALSE(stan::math::is_scal_finite(x));
}

TEST(ErrorHandlingScalar, isScalFinite_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();
  EXPECT_FALSE(stan::math::is_scal_finite(nan));
}
