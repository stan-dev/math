#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, hypergeometric_1f0Double) {
  using stan::math::hypergeometric_1f0;
  using stan::math::inv;

  EXPECT_FLOAT_EQ(4.62962962963, hypergeometric_1f0(3, 0.4));
  EXPECT_FLOAT_EQ(0.510204081633, hypergeometric_1f0(2, -0.4));
  EXPECT_FLOAT_EQ(300.906354890, hypergeometric_1f0(16.0, 0.3));
  EXPECT_FLOAT_EQ(0.531441, hypergeometric_1f0(-6.0, 0.1));
}

TEST(MathFunctions, hypergeometric_1f0_throw) {
  using stan::math::hypergeometric_1f0;

  EXPECT_THROW(hypergeometric_1f0(2.1, 1.0), std::domain_error);
  EXPECT_THROW(hypergeometric_1f0(0.5, 1.5), std::domain_error);
  EXPECT_THROW(hypergeometric_1f0(0.5, -1.0), std::domain_error);
  EXPECT_THROW(hypergeometric_1f0(0.5, -1.5), std::domain_error);
}
