#include <stan/math/prim/scal.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(MathFunctions, round) {
  using stan::math::round;
  EXPECT_FLOAT_EQ(-27, round(-27.3239));
  EXPECT_FLOAT_EQ(-1, round(-0.5));
  EXPECT_FLOAT_EQ(0, round(0));
  EXPECT_FLOAT_EQ(0, round(0.0));
  EXPECT_FLOAT_EQ(1, round(0.5));
  EXPECT_FLOAT_EQ(27, round(27.3239));
}

TEST(MathFunctions, roundNaN) {
  double nan = std::numeric_limits<double>::quiet_NaN();
  EXPECT_PRED1(boost::math::isnan<double>, stan::math::round(nan));
}
