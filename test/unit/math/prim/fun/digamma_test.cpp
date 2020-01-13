#include <stan/math/prim.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, digamma) {
  EXPECT_FLOAT_EQ(boost::math::digamma(0.5), stan::math::digamma(0.5));
  EXPECT_FLOAT_EQ(boost::math::digamma(-1.5), stan::math::digamma(-1.5));
}

TEST(MathFunctions, digamma_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_TRUE(std::isnan(stan::math::digamma(nan)));

  EXPECT_TRUE(std::isnan(stan::math::digamma(-1)));

  EXPECT_TRUE(std::isnormal(stan::math::digamma(1.0E50)));
}
