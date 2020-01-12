#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
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

  EXPECT_TRUE(std::isnan(stan::math::hypot(3.0, nan)));

  EXPECT_TRUE(std::isnan(stan::math::hypot(nan, 3.0)));

  EXPECT_TRUE(std::isnan(stan::math::hypot(nan, nan)));
}
