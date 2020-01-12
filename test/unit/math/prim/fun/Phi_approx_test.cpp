#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, Phi_approx) {
  EXPECT_EQ(0.5, stan::math::Phi_approx(0.0));
  EXPECT_NEAR(stan::math::Phi(0.9), stan::math::Phi_approx(0.9), 0.00014);
  EXPECT_NEAR(stan::math::Phi(-5.0), stan::math::Phi_approx(-5.0), 0.00014);
}

TEST(MathFunctions, Phi_approx_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_TRUE(std::isnan(stan::math::Phi_approx(nan)));
}
