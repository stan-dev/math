#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, log_rising_factorial) {
  using stan::math::log_rising_factorial;

  EXPECT_FLOAT_EQ(std::log(120.0), log_rising_factorial(4.0, 3));
  EXPECT_FLOAT_EQ(std::log(360.0), log_rising_factorial(3.0, 4));
  EXPECT_THROW(log_rising_factorial(-1, 4), std::domain_error);
}

TEST(MathFunctions, log_rising_factorial_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_TRUE(std::isnan(stan::math::log_rising_factorial(nan, 3)));
}
