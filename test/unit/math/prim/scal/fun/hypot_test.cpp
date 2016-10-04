#include <stan/math/prim/scal.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(MathFunctions, hypotDouble) {
  using stan::math::hypot;

  EXPECT_FLOAT_EQ(5.0, hypot(3, 4));
  EXPECT_FLOAT_EQ(5.0, hypot(3, 4.0));
  EXPECT_FLOAT_EQ(5.0, hypot(3.0, 4));
  EXPECT_FLOAT_EQ(5.0, hypot(3.0, 4.0));

  EXPECT_FLOAT_EQ(0.0, hypot(0, 0));
  EXPECT_FLOAT_EQ(0.0, hypot(0, 0.0));
  EXPECT_FLOAT_EQ(0.0, hypot(0.0, 0));
  EXPECT_FLOAT_EQ(0.0, hypot(0.0, 0.0));
}

TEST(MathFunctions, hypotInf) {
  using stan::math::hypot;
  double inf = std::numeric_limits<double>::infinity();

  // promotes results to double
  EXPECT_FLOAT_EQ(inf, hypot(inf, 1));
  EXPECT_FLOAT_EQ(inf, hypot(1, inf));
  EXPECT_FLOAT_EQ(inf, hypot(inf, inf));

  EXPECT_FLOAT_EQ(inf, hypot(-inf, 1));
  EXPECT_FLOAT_EQ(inf, hypot(1, -inf));
  EXPECT_FLOAT_EQ(inf, hypot(-inf, -inf));
  EXPECT_FLOAT_EQ(inf, hypot(-inf, inf));
  EXPECT_FLOAT_EQ(inf, hypot(inf, -inf));
}

TEST(MathFunctions, hypotNaN) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_PRED1(boost::math::isnan<double>,
               stan::math::hypot(3.0, nan));

  EXPECT_PRED1(boost::math::isnan<double>,
               stan::math::hypot(nan, 3.0));

  EXPECT_PRED1(boost::math::isnan<double>,
               stan::math::hypot(nan, nan));
}
