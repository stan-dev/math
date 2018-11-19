#include <stan/math/prim/scal.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <gtest/gtest.h>
#include <limits>

#include <iostream>

// as it just delegates

TEST(MathFunctions, lgamma_finite) {
  EXPECT_FLOAT_EQ(-0.0853740900033158, stan::math::lgamma(1.2));
  EXPECT_FLOAT_EQ(1.57917603403998, stan::math::lgamma(-1.2));
}
TEST(MathFunctions, lgamma_non_pos_int) {
  EXPECT_PRED1(boost::math::isinf<double>, stan::math::lgamma(0));
  EXPECT_PRED1(boost::math::isinf<double>, stan::math::lgamma(0.0));
  EXPECT_PRED1(boost::math::isinf<double>, stan::math::lgamma(-2));
  EXPECT_PRED1(boost::math::isinf<double>, stan::math::lgamma(-2.0));
}

TEST(MathFunctions, lgammaStanMathUsing) { using stan::math::lgamma; }

TEST(MathFunctions, lgamma_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();
  EXPECT_PRED1(boost::math::isnan<double>, stan::math::lgamma(nan));
}
