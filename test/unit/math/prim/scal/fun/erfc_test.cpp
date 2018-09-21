#include <stan/math/prim/scal.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <stdexcept>

TEST(MathFunctions, erfc) {
  using stan::math::erfc;
  EXPECT_FLOAT_EQ(1.983790458590775, erfc(-1.7));
  EXPECT_FLOAT_EQ(1, erfc(0));
  EXPECT_FLOAT_EQ(0.9774354253081551, erfc(0.02));
  EXPECT_FLOAT_EQ(3.558629930076841e-07, erfc(3.6));
}

TEST(MathFunctions, erfcOverflow) {
  EXPECT_FLOAT_EQ(2, erfc(-100));
  EXPECT_FLOAT_EQ(0, erfc(100));
}

TEST(MathFunctions, erfcNan) {
  double nan = std::numeric_limits<double>::quiet_NaN();
  EXPECT_PRED1(boost::math::isnan<double>, stan::math::erfc(nan));
}
