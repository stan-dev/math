#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(ProbGumbel, log_matches_lpdf) {
  double y = 0.8;
  double mu = 1.1;
  double sigma = 2.3;

  EXPECT_FLOAT_EQ((stan::math::gumbel_lpdf(y, mu, sigma)),
                  (stan::math::gumbel_log(y, mu, sigma)));
  EXPECT_FLOAT_EQ((stan::math::gumbel_lpdf<true>(y, mu, sigma)),
                  (stan::math::gumbel_log<true>(y, mu, sigma)));
  EXPECT_FLOAT_EQ((stan::math::gumbel_lpdf<false>(y, mu, sigma)),
                  (stan::math::gumbel_log<false>(y, mu, sigma)));
  EXPECT_FLOAT_EQ((stan::math::gumbel_lpdf<true, double, double, double>(y, mu, sigma)),
                  (stan::math::gumbel_log<true, double, double, double>(y, mu, sigma)));
  EXPECT_FLOAT_EQ((stan::math::gumbel_lpdf<false, double, double, double>(y, mu, sigma)),
                  (stan::math::gumbel_log<false, double, double, double>(y, mu, sigma)));
  EXPECT_FLOAT_EQ((stan::math::gumbel_lpdf<double, double, double>(y, mu, sigma)),
                  (stan::math::gumbel_log<double, double, double>(y, mu, sigma)));
}
