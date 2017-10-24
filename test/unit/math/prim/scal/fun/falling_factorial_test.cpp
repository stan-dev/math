#include <stan/math/prim/scal.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(MathFunctions, falling_factorial) {
  using stan::math::falling_factorial;

  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_FLOAT_EQ(24, falling_factorial(4.0, 3));
  EXPECT_FLOAT_EQ(120, falling_factorial(5.0, 4));
  EXPECT_THROW(falling_factorial(1, -4), std::domain_error);
  EXPECT_THROW(falling_factorial(nan, 1), std::domain_error);
}
