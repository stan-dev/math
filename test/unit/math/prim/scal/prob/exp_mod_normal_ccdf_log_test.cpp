#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(ProbExpModNormal, ccdf_log_matches_lccdf) {
  double y = 0.8;
  double mu = 1.1;
  double sigma = 2.3;
  double lambda = 0.5;
  
  EXPECT_FLOAT_EQ((stan::math::exp_mod_normal_lccdf(y, mu, lambda, sigma)),
                  (stan::math::exp_mod_normal_ccdf_log(y, mu, lambda, sigma)));
  EXPECT_FLOAT_EQ((stan::math::exp_mod_normal_lccdf<double, double, double, double>(y, mu, lambda, sigma)),
                  (stan::math::exp_mod_normal_ccdf_log<double, double, double, double>(y, mu, lambda, sigma)));
}
