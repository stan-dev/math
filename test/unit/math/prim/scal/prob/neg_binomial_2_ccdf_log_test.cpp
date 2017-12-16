#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions.hpp>

TEST(ProbNegBinomial2, ccdf_log_matches_lccdf) {
  double y = 0.8;
  double mu = 1.1;
  double phi = 2.3;

  EXPECT_FLOAT_EQ((stan::math::neg_binomial_2_lccdf(y, mu, phi)),
                  (stan::math::neg_binomial_2_ccdf_log(y, mu, phi)));
  EXPECT_FLOAT_EQ(
      (stan::math::neg_binomial_2_lccdf<double, double, double>(y, mu, phi)),
      (stan::math::neg_binomial_2_ccdf_log<double, double, double>(y, mu,
                                                                   phi)));
}
