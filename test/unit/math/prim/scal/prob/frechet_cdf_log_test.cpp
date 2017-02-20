#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(ProbFrechet, cdf_log_matches_lcdf) {
  double y = 0.3;
  double alpha = 2;
  double sigma = 1.5;

  EXPECT_FLOAT_EQ((stan::math::frechet_lcdf(y, alpha, sigma)),
                  (stan::math::frechet_cdf_log(y, alpha, sigma)));
  EXPECT_FLOAT_EQ((stan::math::frechet_lcdf<double, double>(y, alpha, sigma)),
                  (stan::math::frechet_cdf_log<double, double>(y, alpha, sigma)));
}
