#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(ProbPoisson, cdf_log_matches_lcdf) {
  int y = 3;
  double lambda = 2.3;

  EXPECT_FLOAT_EQ((stan::math::poisson_lcdf(y, lambda)),
                  (stan::math::poisson_cdf_log(y, lambda)));
  EXPECT_FLOAT_EQ((stan::math::poisson_lcdf<int, double>(y, lambda)),
                  (stan::math::poisson_cdf_log<int, double>(y, lambda)));
}
