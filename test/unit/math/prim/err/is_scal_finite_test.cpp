#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(ErrorHandlingScalar, isScalFinite) {
  using stan::math::is_scal_finite;
  double x = 0;
  EXPECT_TRUE(is_scal_finite(x));

  x = std::numeric_limits<double>::infinity();
  EXPECT_FALSE(is_scal_finite(x));

  x = -std::numeric_limits<double>::infinity();
  EXPECT_FALSE(is_scal_finite(x));

  x = std::numeric_limits<double>::quiet_NaN();
  EXPECT_FALSE(is_scal_finite(x));
}

TEST(ErrorHandlingScalar, isScalFinite_nan) {
  using stan::math::is_scal_finite;
  double nan = std::numeric_limits<double>::quiet_NaN();
  EXPECT_FALSE(is_scal_finite(nan));
}
