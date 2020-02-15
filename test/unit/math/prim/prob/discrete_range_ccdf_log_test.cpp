#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(ProbDiscreteRange, cdf_log_matches_lccdf) {
  using stan::math::discrete_range_ccdf_log;
  using stan::math::discrete_range_lccdf;
  int y = 1;
  int lower = 0;
  int upper = 9;

  EXPECT_FLOAT_EQ((discrete_range_lccdf(y, lower, upper)),
                  (discrete_range_ccdf_log(y, lower, upper)));
  EXPECT_FLOAT_EQ((discrete_range_lccdf<double>(y, lower, upper)),
                  (discrete_range_ccdf_log<double>(y, lower, upper)));
}
