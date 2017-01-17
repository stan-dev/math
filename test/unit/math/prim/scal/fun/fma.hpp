#include <stan/math/prim/scal.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <gtest/gtest.h>

TEST(MathFunctions, fma_double) {
  using stan::math::fma;
  EXPECT_FLOAT_EQ(1.0, fma(3.0, 2.0, -5));
  EXPECT_FLOAT_EQ(0.0, fma(2.0, 3.0, -6));
  EXPECT_FLOAT_EQ(46.9, fma(4.5, 2.0, 37.9));
}

TEST(MathFunctions, fma_int) {
  using stan::math::fma;
  EXPECT_FLOAT_EQ(11.0, fma(int(3),int(2),int(5)));
}

TEST(MathFunctions, fma_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();
  EXPECT_PRED1(boost::math::isnan<double>,
               stan::math::fma(nan, 3.0, 2.7));
  EXPECT_PRED1(boost::math::isnan<double>,
               stan::math::fma(3.0, nan, 1.5));
  EXPECT_PRED1(boost::math::isnan<double>,
               stan::math::fma(2, -8.2, nan));
}
