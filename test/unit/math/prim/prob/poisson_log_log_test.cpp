#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(ProbPoissonLog, log_matches_lpmf) {
  int y = 3;
  double alpha = -0.3;

  EXPECT_FLOAT_EQ((stan::math::poisson_log_lpmf(y, alpha)),
                  (stan::math::poisson_log_log(y, alpha)));
  EXPECT_FLOAT_EQ((stan::math::poisson_log_lpmf<true>(y, alpha)),
                  (stan::math::poisson_log_log<true>(y, alpha)));
  EXPECT_FLOAT_EQ((stan::math::poisson_log_lpmf<false>(y, alpha)),
                  (stan::math::poisson_log_log<false>(y, alpha)));
  EXPECT_FLOAT_EQ((stan::math::poisson_log_lpmf<true, int, double>(y, alpha)),
                  (stan::math::poisson_log_log<true, int, double>(y, alpha)));
  EXPECT_FLOAT_EQ((stan::math::poisson_log_lpmf<false, int, double>(y, alpha)),
                  (stan::math::poisson_log_log<false, int, double>(y, alpha)));
  EXPECT_FLOAT_EQ((stan::math::poisson_log_lpmf<int, double>(y, alpha)),
                  (stan::math::poisson_log_log<int, double>(y, alpha)));
}
