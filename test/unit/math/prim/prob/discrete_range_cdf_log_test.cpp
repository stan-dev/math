#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(ProbDiscreteRange, cdf_log_matches_lcdf) {
  using stan::math::discrete_range_cdf_log;
  using stan::math::discrete_range_lcdf;
  int y = 1;
  int lower = 0;
  int upper = 9;

  EXPECT_FLOAT_EQ((discrete_range_lcdf(y, lower, upper)),
                  (discrete_range_cdf_log(y, lower, upper)));
  EXPECT_FLOAT_EQ((discrete_range_lcdf<double>(y, lower, upper)),
                  (discrete_range_cdf_log<double>(y, lower, upper)));
}
