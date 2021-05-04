#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(ProbDiscreteRange, log_matches_lpmf) {
  int y = 3;
  int lower = -2;
  int upper = 5;

  EXPECT_FLOAT_EQ((stan::math::discrete_range_lpmf(y, lower, upper)),
                  (stan::math::discrete_range_log(y, lower, upper)));
  EXPECT_FLOAT_EQ((stan::math::discrete_range_lpmf<true>(y, lower, upper)),
                  (stan::math::discrete_range_log<true>(y, lower, upper)));
  EXPECT_FLOAT_EQ((stan::math::discrete_range_lpmf<false>(y, lower, upper)),
                  (stan::math::discrete_range_log<false>(y, lower, upper)));
  EXPECT_FLOAT_EQ(
      (stan::math::discrete_range_lpmf<true, int, int, int>(y, lower, upper)),
      (stan::math::discrete_range_log<true, int, int, int>(y, lower, upper)));
  EXPECT_FLOAT_EQ(
      (stan::math::discrete_range_lpmf<false, int, int, int>(y, lower, upper)),
      (stan::math::discrete_range_log<false, int, int, int>(y, lower, upper)));
  EXPECT_FLOAT_EQ(
      (stan::math::discrete_range_lpmf<int, int, int>(y, lower, upper)),
      (stan::math::discrete_range_log<int, int, int>(y, lower, upper)));
}
