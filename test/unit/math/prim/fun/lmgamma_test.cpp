#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, lmgamma) {
  unsigned int k = 1;
  double x = 2.5;
  double result = k * (k - 1) * log(boost::math::constants::pi<double>()) / 4.0;
  // j = 1
  result += lgamma(x);
  EXPECT_FLOAT_EQ(result, stan::math::lmgamma(k, x));

  k = 2;
  x = 3.0;
  result = k * (k - 1) * log(boost::math::constants::pi<double>()) / 4.0;
  // j = 1
  result += lgamma(x);
  // j = 2
  result += lgamma(x + (1.0 - 2.0) / 2.0);
  EXPECT_FLOAT_EQ(result, stan::math::lmgamma(k, x));
}

TEST(MathFunctions, lmgamma_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_TRUE(std::isnan(stan::math::lmgamma(2, nan)));
}
