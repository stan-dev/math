#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(ProbCauchy, log_matches_lpdf) {
  double y = 2;
  double mu = 1.5;
  double sigma = 0.5;

  EXPECT_FLOAT_EQ((stan::math::cauchy_lpdf(y, mu, sigma)),
                  (stan::math::cauchy_log(y, mu, sigma)));
  EXPECT_FLOAT_EQ((stan::math::cauchy_lpdf<true>(y, mu, sigma)),
                  (stan::math::cauchy_log<true>(y, mu, sigma)));
  EXPECT_FLOAT_EQ((stan::math::cauchy_lpdf<false>(y, mu, sigma)),
                  (stan::math::cauchy_log<false>(y, mu, sigma)));
  EXPECT_FLOAT_EQ((stan::math::cauchy_lpdf<true, double, double, double>(y, mu, sigma)),
                  (stan::math::cauchy_log<true, double, double, double>(y, mu, sigma)));
  EXPECT_FLOAT_EQ((stan::math::cauchy_lpdf<false, double, double, double>(y, mu, sigma)),
                  (stan::math::cauchy_log<false, double, double, double>(y, mu, sigma)));
  EXPECT_FLOAT_EQ((stan::math::cauchy_lpdf<double, double, double>(y, mu, sigma)),
                  (stan::math::cauchy_log<double, double, double>(y, mu, sigma)));
}
