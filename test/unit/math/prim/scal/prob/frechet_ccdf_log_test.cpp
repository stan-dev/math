#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(ProbFrechet, ccdf_log_matches_lccdf) {
  double y = 0.3;
  double alpha = 2;
  double sigma = 1.5;

  EXPECT_FLOAT_EQ((stan::math::frechet_lccdf(y, alpha, sigma)),
                  (stan::math::frechet_ccdf_log(y, alpha, sigma)));
  EXPECT_FLOAT_EQ((stan::math::frechet_lccdf<double, double>(y, alpha, sigma)),
                  (stan::math::frechet_ccdf_log<double, double>(y, alpha, sigma)));
}
