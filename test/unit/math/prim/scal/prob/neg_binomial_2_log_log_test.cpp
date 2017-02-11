#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(ProbNegBinomial2Log, log_matches_lpmf) {
  double y = 0.8;
  double eta = 1.1;
  double phi = 2.3;

  EXPECT_FLOAT_EQ((stan::math::neg_binomial_2_log_lpmf(y, eta, phi)),
                  (stan::math::neg_binomial_2_log_log(y, eta, phi)));
  EXPECT_FLOAT_EQ((stan::math::neg_binomial_2_log_lpmf<true>(y, eta, phi)),
                  (stan::math::neg_binomial_2_log_log<true>(y, eta, phi)));
  EXPECT_FLOAT_EQ((stan::math::neg_binomial_2_log_lpmf<false>(y, eta, phi)),
                  (stan::math::neg_binomial_2_log_log<false>(y, eta, phi)));
  EXPECT_FLOAT_EQ((stan::math::neg_binomial_2_log_lpmf<true, double, double, double>(y, eta, phi)),
                  (stan::math::neg_binomial_2_log_log<true, double, double, double>(y, eta, phi)));
  EXPECT_FLOAT_EQ((stan::math::neg_binomial_2_log_lpmf<false, double, double, double>(y, eta, phi)),
                  (stan::math::neg_binomial_2_log_log<false, double, double, double>(y, eta, phi)));
  EXPECT_FLOAT_EQ((stan::math::neg_binomial_2_log_lpmf<double, double, double>(y, eta, phi)),
                  (stan::math::neg_binomial_2_log_log<double, double, double>(y, eta, phi)));
}
