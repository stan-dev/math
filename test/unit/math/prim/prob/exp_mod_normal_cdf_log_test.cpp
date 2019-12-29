#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(ProbFoo, cdf_log_matches_lcdf) {
  double y = 0.8;
  double mu = 1.1;
  double sigma = 2.3;
  double lambda = 0.5;

  EXPECT_FLOAT_EQ((stan::math::exp_mod_normal_lcdf(y, mu, lambda, sigma)),
                  (stan::math::exp_mod_normal_cdf_log(y, mu, lambda, sigma)));
  EXPECT_FLOAT_EQ(
      (stan::math::exp_mod_normal_lcdf<double, double, double, double>(
          y, mu, lambda, sigma)),
      (stan::math::exp_mod_normal_cdf_log<double, double, double, double>(
          y, mu, lambda, sigma)));
}
