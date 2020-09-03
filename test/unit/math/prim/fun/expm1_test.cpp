#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <stdexcept>

TEST(MathFunctions, expm1) {
  using stan::math::expm1;
  using std::exp;
  EXPECT_FLOAT_EQ(exp(-14.2) - 1, expm1(-14.2));
  EXPECT_FLOAT_EQ(exp(0) - 1, expm1(0));
  EXPECT_FLOAT_EQ(exp(172.987) - 1, expm1(172.987));
  EXPECT_FLOAT_EQ(-1, expm1(-std::numeric_limits<double>::infinity()));
}

TEST(MathFunctions, expm1_inf_return) {
  EXPECT_EQ(std::numeric_limits<double>::infinity(),
            stan::math::expm1(std::numeric_limits<double>::infinity()));
}

TEST(MathFunctions, expm1_nan) {
  using stan::math::expm1;
  double nan = std::numeric_limits<double>::quiet_NaN();
  EXPECT_TRUE(std::isnan(stan::math::expm1(nan)));
}
