#include <stan/math/prim/scal.hpp>
#include <limits>
#include <stdexcept>
#include <cmath>
#include <boost/math/special_functions/fpclassify.hpp>
#include <gtest/gtest.h>

TEST(MathFunctions, acosh) {
  using stan::math::acosh;
  EXPECT_FLOAT_EQ(0, acosh(1));
  EXPECT_FLOAT_EQ(0.96242365, acosh(1.5));
  EXPECT_FLOAT_EQ(3.0797991, acosh(10.9));
}

TEST(MathFunctions, acosh_exception) {
  using stan::math::acosh;
  EXPECT_THROW(acosh(0.5), std::domain_error);
  EXPECT_THROW(acosh(std::numeric_limits<double>::infinity()), std::overflow_error);
}

TEST(MathFunctions, log1p_nan) {
  using stan::math::acosh;
  EXPECT_PRED1(boost::math::isnan<double>,
               stan::math::acosh(std::numeric_limits<double>::quiet_NaN()));
}
