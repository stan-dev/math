#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(ProbNegBinomial, ccdf_log_matches_lccdf) {
  int y = 3;
  double alpha = 1.1;
  double beta = 2.3;

  EXPECT_FLOAT_EQ((stan::math::neg_binomial_lccdf(y, alpha, beta)),
                  (stan::math::neg_binomial_ccdf_log(y, alpha, beta)));
  EXPECT_FLOAT_EQ(
      (stan::math::neg_binomial_lccdf<int, double, double>(y, alpha, beta)),
      (stan::math::neg_binomial_ccdf_log<int, double, double>(y, alpha, beta)));
}
