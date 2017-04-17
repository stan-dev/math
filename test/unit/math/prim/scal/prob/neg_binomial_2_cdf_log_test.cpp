#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

 TEST(ProbNegBinomial2, cdf_log_matches_lcdf) {
  double y = 0.8;
  double mu = 1.1;
  double phi = 2.3;

  EXPECT_FLOAT_EQ((stan::math::neg_binomial_2_lcdf(y, mu, phi)),
                  (stan::math::neg_binomial_2_cdf_log(y, mu, phi)));
  EXPECT_FLOAT_EQ((stan::math::neg_binomial_2_lcdf<double, double, double>(y, mu, phi)),
                  (stan::math::neg_binomial_2_cdf_log<double, double, double>(y, mu, phi)));
}
