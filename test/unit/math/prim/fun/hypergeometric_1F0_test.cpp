#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, hypergeometric_1f0Double) {
  using stan::math::hypergeometric_1f0;
  using stan::math::inv;

  EXPECT_FLOAT_EQ(-inv(27.0), hypergeometric_1f0(3, 4));
  EXPECT_FLOAT_EQ(inv(25.0), hypergeometric_1f0(2, -4.0));
  EXPECT_FLOAT_EQ(inv(65536), hypergeometric_1f0(16.0, 3));
  EXPECT_FLOAT_EQ(531441.0, hypergeometric_1f0(-6.0, 10.0));
}

TEST(MathFunctions, hypergeometric_1f0_throw) {
  using stan::math::hypergeometric_1f0;

  EXPECT_THROW(hypergeometric_1f0(2.1, 1.0), std::domain_error);
  EXPECT_THROW(hypergeometric_1f0(0.5, 1.5), std::domain_error);
}
