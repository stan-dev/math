#include <stan/math/prim/scal.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(MathFunctions, rising_factorial) {
  using stan::math::rising_factorial;

  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_FLOAT_EQ(120, rising_factorial(4.0, 3));
  EXPECT_FLOAT_EQ(360, rising_factorial(3.0, 4));
  EXPECT_THROW(rising_factorial(1, -4), std::domain_error);
  EXPECT_THROW(rising_factorial(nan, 1), std::domain_error);
}
