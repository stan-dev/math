#include <stan/math/prim/scal.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(MathFunctions, exp2_double) {
  using stan::math::exp2;

  EXPECT_FLOAT_EQ(1.0, exp2(0.0));
  EXPECT_FLOAT_EQ(2.0, exp2(1.0));
  EXPECT_FLOAT_EQ(4.0, exp2(2.0));
  EXPECT_FLOAT_EQ(8.0, exp2(3.0));

  EXPECT_FLOAT_EQ(0.5, exp2(-1.0));
  EXPECT_FLOAT_EQ(0.25, exp2(-2.0));
  EXPECT_FLOAT_EQ(0.125, exp2(-3.0));
}

TEST(MathFunctions, exp2_int) {
  using stan::math::exp2;

  // promotes results to double
  EXPECT_FLOAT_EQ(1.0, exp2(static_cast<int>(0)));
  EXPECT_FLOAT_EQ(2.0, exp2(static_cast<int>(1)));
  EXPECT_FLOAT_EQ(4.0, exp2(static_cast<int>(2)));
  EXPECT_FLOAT_EQ(8.0, exp2(static_cast<int>(3)));

  EXPECT_FLOAT_EQ(0.5, exp2(static_cast<int>(-1)));
  EXPECT_FLOAT_EQ(0.25, exp2(static_cast<int>(-2)));
  EXPECT_FLOAT_EQ(0.125, exp2(static_cast<int>(-3)));
}

TEST(MathFunctions, exp2_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_PRED1(boost::math::isnan<double>, stan::math::exp2(nan));
}
