#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(ProbNegBinomial, log_matches_lpmf) {
  int y = 3;
  double alpha = 1.1;
  double beta = 2.3;

  EXPECT_FLOAT_EQ((stan::math::neg_binomial_lpmf(y, alpha, beta)),
                  (stan::math::neg_binomial_log(y, alpha, beta)));
  EXPECT_FLOAT_EQ((stan::math::neg_binomial_lpmf<true>(y, alpha, beta)),
                  (stan::math::neg_binomial_log<true>(y, alpha, beta)));
  EXPECT_FLOAT_EQ((stan::math::neg_binomial_lpmf<false>(y, alpha, beta)),
                  (stan::math::neg_binomial_log<false>(y, alpha, beta)));
  EXPECT_FLOAT_EQ((stan::math::neg_binomial_lpmf<true, int, double, double>(y, alpha, beta)),
                  (stan::math::neg_binomial_log<true, int, double, double>(y, alpha, beta)));
  EXPECT_FLOAT_EQ((stan::math::neg_binomial_lpmf<false, int, double, double>(y, alpha, beta)),
                  (stan::math::neg_binomial_log<false, int, double, double>(y, alpha, beta)));
  EXPECT_FLOAT_EQ((stan::math::neg_binomial_lpmf<int, double, double>(y, alpha, beta)),
                  (stan::math::neg_binomial_log<int, double, double>(y, alpha, beta)));
}
