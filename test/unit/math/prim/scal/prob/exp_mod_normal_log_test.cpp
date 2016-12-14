#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(ProbExpModNormal, log_matches_lpmf) {
  double y = 0.8;
  double mu = 1.1;
  double sigma = 2.3;
  double lambda = 0.5;
  
  EXPECT_FLOAT_EQ((stan::math::exp_mod_normal_lpdf(y, mu, lambda, sigma)),
                  (stan::math::exp_mod_normal_log(y, mu, lambda, sigma)));
  EXPECT_FLOAT_EQ((stan::math::exp_mod_normal_lpdf<true>(y, mu, lambda, sigma)),
                  (stan::math::exp_mod_normal_log<true>(y, mu, lambda, sigma)));
  EXPECT_FLOAT_EQ((stan::math::exp_mod_normal_lpdf<false>(y, mu, lambda, sigma)),
                  (stan::math::exp_mod_normal_log<false>(y, mu, lambda, sigma)));
  EXPECT_FLOAT_EQ((stan::math::exp_mod_normal_lpdf<true, double, double, double, double>(y, mu, lambda, sigma)),
                  (stan::math::exp_mod_normal_log<true, double, double, double, double>(y, mu, lambda, sigma)));
  EXPECT_FLOAT_EQ((stan::math::exp_mod_normal_lpdf<false, double, double, double, double>(y, mu, lambda, sigma)),
                  (stan::math::exp_mod_normal_log<false, double, double, double, double>(y, mu, lambda, sigma)));
  EXPECT_FLOAT_EQ((stan::math::exp_mod_normal_lpdf<double, double, double, double>(y, mu, lambda, sigma)),
                  (stan::math::exp_mod_normal_log<double, double, double, double>(y, mu, lambda, sigma)));
}
