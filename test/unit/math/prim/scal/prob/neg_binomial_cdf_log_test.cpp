#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(ProbNegBinomial, cdf_log_matches_lcdf) {
  int y = 3;
  double alpha = 1.1;
  double beta = 2.3;

  EXPECT_FLOAT_EQ((stan::math::neg_binomial_lcdf(y, alpha, beta)),
                  (stan::math::neg_binomial_cdf_log(y, alpha, beta)));
  EXPECT_FLOAT_EQ(
      (stan::math::neg_binomial_lcdf<int, double, double>(y, alpha, beta)),
      (stan::math::neg_binomial_cdf_log<int, double, double>(y, alpha, beta)));
}
