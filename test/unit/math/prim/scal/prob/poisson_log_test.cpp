#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(ProbPoisson, log_matches_lpmf) {
  int y = 3;
  double lambda = 2.3;

  EXPECT_FLOAT_EQ((stan::math::poisson_lpmf(y, lambda)),
                  (stan::math::poisson_log(y, lambda)));
  EXPECT_FLOAT_EQ((stan::math::poisson_lpmf<true>(y, lambda)),
                  (stan::math::poisson_log<true>(y, lambda)));
  EXPECT_FLOAT_EQ((stan::math::poisson_lpmf<false>(y, lambda)),
                  (stan::math::poisson_log<false>(y, lambda)));
  EXPECT_FLOAT_EQ((stan::math::poisson_lpmf<true>(y, lambda)),
                  (stan::math::poisson_log<true>(y, lambda)));
  EXPECT_FLOAT_EQ((stan::math::poisson_lpmf<false>(y, lambda)),
                  (stan::math::poisson_log<false>(y, lambda)));
  EXPECT_FLOAT_EQ((stan::math::poisson_lpmf(y, lambda)),
                  (stan::math::poisson_log(y, lambda)));
}
