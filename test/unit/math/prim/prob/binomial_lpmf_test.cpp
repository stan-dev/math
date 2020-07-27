#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(ProbBinomial, edge_cases) {
  EXPECT_FLOAT_EQ(stan::math::binomial_lpmf(2, 2, 1.0), 0.0);
  EXPECT_TRUE(stan::math::is_inf(stan::math::binomial_lpmf(0, 2, 1.0)));
  EXPECT_TRUE(stan::math::is_inf(stan::math::binomial_lpmf(2, 2, 0.0)));
  EXPECT_FLOAT_EQ(stan::math::binomial_lpmf(0, 2, 0.0), 0.0);
}
