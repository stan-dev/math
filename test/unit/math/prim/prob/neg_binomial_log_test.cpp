#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(ProbNegBinomial, log_matches_lpmf) {
  int y = 3;
  double alpha = 1.1;
  double beta = 2.3;

  long double alpha2 = 1e11;
  long double beta2 = 1e10;

  EXPECT_FLOAT_EQ(stan::math::neg_binomial_lpmf(0, 1e11, 1e10),
                  -alpha2 * std::log1p(1.0 / beta2));
  EXPECT_FLOAT_EQ(stan::math::neg_binomial_lpmf(1, 1e11, 1e10),
                  std::log(alpha2 - 1.0) - alpha2 * std::log1p(1.0 / beta2)
                      - std::log1p(beta2));
  EXPECT_FLOAT_EQ((stan::math::neg_binomial_lpmf(13, 1e11, 1e10)),
                  -2.6185576442208003);

  EXPECT_FLOAT_EQ((stan::math::neg_binomial_lpmf(y, alpha, beta)),
                  (stan::math::neg_binomial_log(y, alpha, beta)));
  EXPECT_FLOAT_EQ((stan::math::neg_binomial_lpmf<true>(y, alpha, beta)),
                  (stan::math::neg_binomial_log<true>(y, alpha, beta)));
  EXPECT_FLOAT_EQ((stan::math::neg_binomial_lpmf<false>(y, alpha, beta)),
                  (stan::math::neg_binomial_log<false>(y, alpha, beta)));
  EXPECT_FLOAT_EQ((stan::math::neg_binomial_lpmf<true, int, double, double>(
                      y, alpha, beta)),
                  (stan::math::neg_binomial_log<true, int, double, double>(
                      y, alpha, beta)));
  EXPECT_FLOAT_EQ((stan::math::neg_binomial_lpmf<false, int, double, double>(
                      y, alpha, beta)),
                  (stan::math::neg_binomial_log<false, int, double, double>(
                      y, alpha, beta)));
  EXPECT_FLOAT_EQ(
      (stan::math::neg_binomial_lpmf<int, double, double>(y, alpha, beta)),
      (stan::math::neg_binomial_log<int, double, double>(y, alpha, beta)));
}
