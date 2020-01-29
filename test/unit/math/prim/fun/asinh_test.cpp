#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <stdexcept>

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
  double nan = std::numeric_limits<double>::quiet_NaN();
  EXPECT_TRUE(std::isnan(stan::math::asinh(nan)));
}
