#include <stan/math/prim.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, tgamma) {
  using stan::math::tgamma;
  EXPECT_FLOAT_EQ(1.772453850905516, tgamma(0.5));
  EXPECT_FLOAT_EQ(1, tgamma(1));
  EXPECT_FLOAT_EQ(2.423965479935368, tgamma(3.2));
  EXPECT_FLOAT_EQ(6402373705728000, tgamma(19));
}

TEST(MathFunctions, tgammaStanMath) {
  EXPECT_THROW(stan::math::tgamma(0.0), std::domain_error);
}

TEST(MathFunctions, tgamma_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();
  EXPECT_TRUE(std::isnan(stan::math::tgamma(nan)));
}
