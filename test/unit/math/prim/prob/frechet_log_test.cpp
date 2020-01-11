#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(ProbFrechet, log_matches_lpdf) {
  double y = 0.3;
  double alpha = 2;
  double sigma = 1.5;

  EXPECT_FLOAT_EQ((stan::math::frechet_lpdf(y, alpha, sigma)),
                  (stan::math::frechet_log(y, alpha, sigma)));
  EXPECT_FLOAT_EQ((stan::math::frechet_lpdf<true>(y, alpha, sigma)),
                  (stan::math::frechet_log<true>(y, alpha, sigma)));
  EXPECT_FLOAT_EQ((stan::math::frechet_lpdf<false>(y, alpha, sigma)),
                  (stan::math::frechet_log<false>(y, alpha, sigma)));
  EXPECT_FLOAT_EQ(
      (stan::math::frechet_lpdf<true, double, double>(y, alpha, sigma)),
      (stan::math::frechet_log<true, double, double>(y, alpha, sigma)));
  EXPECT_FLOAT_EQ(
      (stan::math::frechet_lpdf<false, double, double>(y, alpha, sigma)),
      (stan::math::frechet_log<false, double, double>(y, alpha, sigma)));
  EXPECT_FLOAT_EQ((stan::math::frechet_lpdf<double, double>(y, alpha, sigma)),
                  (stan::math::frechet_log<double, double>(y, alpha, sigma)));
}
