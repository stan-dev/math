#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, lbeta) {
  using stan::math::lbeta;

  EXPECT_FLOAT_EQ(0.0, lbeta(1.0, 1.0));
  EXPECT_FLOAT_EQ(2.981361, lbeta(0.1, 0.1));
  EXPECT_FLOAT_EQ(-4.094345, lbeta(3.0, 4.0));
  EXPECT_FLOAT_EQ(-4.094345, lbeta(4.0, 3.0));
}

TEST(MathFunctions, lbeta_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_TRUE(std::isnan(stan::math::lbeta(nan, 1.0)));

  EXPECT_TRUE(std::isnan(stan::math::lbeta(1.0, nan)));

  EXPECT_TRUE(std::isnan(stan::math::lbeta(nan, nan)));
}
