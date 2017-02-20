#include <stan/math/prim/scal.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <gtest/gtest.h>

TEST(MathFunctions, trunc) {
  using stan::math::trunc;
  EXPECT_FLOAT_EQ(-27, trunc(-27.8239));
  EXPECT_FLOAT_EQ(-1, trunc(-1.5));
  EXPECT_FLOAT_EQ(0, trunc(0));
  EXPECT_FLOAT_EQ(0, trunc(0.0));
  EXPECT_FLOAT_EQ(0, trunc(0.5));
  EXPECT_FLOAT_EQ(27, trunc(27.3239));
} 

TEST(MathFunctions, truncNaN) {
  double nan = std::numeric_limits<double>::quiet_NaN();
  EXPECT_PRED1(boost::math::isnan<double>,
               stan::math::trunc(nan));
}
