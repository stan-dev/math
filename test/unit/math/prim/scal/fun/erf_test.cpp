#include <stan/math/prim/scal.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <stdexcept>

TEST(MathFunctions, erf) {
  using stan::math::erf;
  EXPECT_FLOAT_EQ(-0.3286267594591274, erf(-0.3));
  EXPECT_FLOAT_EQ(0, erf(0));
  EXPECT_FLOAT_EQ(0.9999939742388482, erf(3.2));
}

TEST(MathFunctions, erfOverflow) {
  EXPECT_FLOAT_EQ(-1, erf(-100));
  EXPECT_FLOAT_EQ(1, erf(100));
}

TEST(MathFunctions, erfNan) {
  double nan = std::numeric_limits<double>::quiet_NaN();
  EXPECT_PRED1(boost::math::isnan<double>, stan::math::erf(nan));
}
