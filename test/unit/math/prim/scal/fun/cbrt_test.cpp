#include <stan/math/prim/scal.hpp>
#include <limits>
#include <stdexcept>
#include <cmath>
#include <boost/math/special_functions/fpclassify.hpp>
#include <gtest/gtest.h>

TEST(MathFunctions, cbrt) {
  using stan::math::cbrt;
  EXPECT_FLOAT_EQ(-2.0, cbrt(-8.0));
  EXPECT_FLOAT_EQ(-1.392476650083834, cbrt(-2.7));
  EXPECT_FLOAT_EQ(0, cbrt(0));
  EXPECT_FLOAT_EQ(2.0, cbrt(8.0));

}

TEST(MathFunctions, cbrt_inf_return) {
  EXPECT_EQ(-std::numeric_limits<double>::infinity(),
            stan::math::cbrt(-std::numeric_limits<double>::infinity()));
  EXPECT_EQ(std::numeric_limits<double>::infinity(),
            stan::math::cbrt(std::numeric_limits<double>::infinity()));
}

TEST(MathFunctions, cbrt_nan) {
  using stan::math::cbrt;
  EXPECT_PRED1(boost::math::isnan<double>,
               stan::math::cbrt(std::numeric_limits<double>::quiet_NaN()));
}
