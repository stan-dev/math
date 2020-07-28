#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(ProbBinomial, edge_cases_zero) {
  EXPECT_FLOAT_EQ(stan::math::binomial_lpmf(0, 0, 0.0), 0.0);
  EXPECT_FLOAT_EQ(stan::math::binomial_lpmf(0, 0, 0.5), 0.0);
  EXPECT_FLOAT_EQ(stan::math::binomial_lpmf(0, 0, 1.0), 0.0);
}

TEST(ProbBinomial, edge_cases_one) {
  EXPECT_FLOAT_EQ(stan::math::binomial_lpmf(0, 1, 0.0), 0.0);
  EXPECT_TRUE(stan::math::is_inf(stan::math::binomial_lpmf(1, 1, 0.0)));
  EXPECT_FLOAT_EQ(stan::math::binomial_lpmf(0, 1, 0.5), std::log(0.5));
  EXPECT_FLOAT_EQ(stan::math::binomial_lpmf(1, 1, 0.5), std::log(0.5));
  EXPECT_TRUE(stan::math::is_inf(stan::math::binomial_lpmf(0, 1, 1.0)));
  EXPECT_FLOAT_EQ(stan::math::binomial_lpmf(1, 1, 1.0), 0.0);
}

TEST(ProbBinomial, edge_cases_two) {
  EXPECT_FLOAT_EQ(stan::math::binomial_lpmf(0, 2, 0.0), 0.0);
  EXPECT_TRUE(stan::math::is_inf(stan::math::binomial_lpmf(1, 2, 0.0)));
  EXPECT_TRUE(stan::math::is_inf(stan::math::binomial_lpmf(2, 2, 0.0)));
  EXPECT_FLOAT_EQ(stan::math::binomial_lpmf(0, 2, 0.5), std::log(0.25));
  EXPECT_FLOAT_EQ(stan::math::binomial_lpmf(1, 2, 0.5), std::log(0.5));
  EXPECT_FLOAT_EQ(stan::math::binomial_lpmf(2, 2, 0.5), std::log(0.25));
  EXPECT_TRUE(stan::math::is_inf(stan::math::binomial_lpmf(0, 2, 1.0)));
  EXPECT_TRUE(stan::math::is_inf(stan::math::binomial_lpmf(1, 2, 1.0)));
  EXPECT_FLOAT_EQ(stan::math::binomial_lpmf(2, 2, 1.0), 0.0);
}
