#include <stan/math/prim/scal.hpp>
#include <limits>
#include <stdexcept>
#include <cmath>
#include <boost/math/special_functions/fpclassify.hpp>
#include <gtest/gtest.h>

TEST(MathFunctions, asinh) {
  using stan::math::asinh;
  EXPECT_FLOAT_EQ(-3.347626679085641, asinh(-14.2));
  EXPECT_FLOAT_EQ(0, asinh(0));
  EXPECT_FLOAT_EQ(5.846371981978736, asinh(172.987));
}

TEST(MathFunctions, asinh_inf_return) {
  EXPECT_EQ(-std::numeric_limits<double>::infinity(),
            stan::math::asinh(-std::numeric_limits<double>::infinity()));
  EXPECT_EQ(std::numeric_limits<double>::infinity(),
            stan::math::asinh(std::numeric_limits<double>::infinity()));
}

TEST(MathFunctions, asinh_nan) {
  using stan::math::asinh;
  EXPECT_PRED1(boost::math::isnan<double>,
               stan::math::asinh(std::numeric_limits<double>::quiet_NaN()));
}
