#include <stan/math/prim/scal.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(MathFunctions, log_falling_factorial) {
  using stan::math::log_falling_factorial;

  EXPECT_FLOAT_EQ(std::log(24), log_falling_factorial(4.0, 3));
  EXPECT_FLOAT_EQ(std::log(120), log_falling_factorial(5.0, 4));
  EXPECT_FLOAT_EQ(std::log(144.8261), log_falling_factorial(6.1, 3.1));
  EXPECT_THROW(log_falling_factorial(-1, 4), std::domain_error);
}

TEST(MathFunctions, log_falling_factorial_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_PRED1(boost::math::isnan<double>,
               stan::math::log_falling_factorial(nan, 3));
}
