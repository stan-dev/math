#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(ProbPoisson, ccdf_log_matches_lccdf) {
  int y = 3;
  double lambda = 2.3;

  EXPECT_FLOAT_EQ((stan::math::poisson_lccdf(y, lambda)),
                  (stan::math::poisson_ccdf_log(y, lambda)));
  EXPECT_FLOAT_EQ((stan::math::poisson_lccdf<int, double>(y, lambda)),
                  (stan::math::poisson_ccdf_log<int, double>(y, lambda)));
}
