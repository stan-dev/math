#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(ProbDiscreteRange, cdf_log_matches_lcdf) {
  using stan::math::discrete_range_cdf_log;
  using stan::math::discrete_range_lcdf;

  for (int lower = 0; lower < 5; ++lower) {
    for (int upper = lower; upper < 5; ++upper) {
      for (int y = lower; y <= upper; ++y) {
        EXPECT_FLOAT_EQ((discrete_range_lcdf(y, lower, upper)),
                        (discrete_range_cdf_log(y, lower, upper)));
        EXPECT_FLOAT_EQ((discrete_range_lcdf<double>(y, lower, upper)),
                        (discrete_range_cdf_log<double>(y, lower, upper)));
      }
    }
  }
}

TEST(ProbDiscreteRange, lcdf_boundaries) {
  using stan::math::discrete_range_lcdf;
  int lower = 1;
  int upper = 5;

  EXPECT_FLOAT_EQ(std::log(1.0 / 5), discrete_range_lcdf(lower, lower, upper));
  EXPECT_FLOAT_EQ(0, discrete_range_lcdf(upper, lower, upper));
}

TEST(ProbDiscreteRange, lcdf_out_of_support) {
  using stan::math::discrete_range_lcdf;
  int lower = 1;
  int upper = 5;

  EXPECT_FLOAT_EQ(stan::math::LOG_ZERO,
                  discrete_range_lcdf(lower - 1, lower, upper));
  EXPECT_FLOAT_EQ(0.0, discrete_range_lcdf(upper + 1, lower, upper));
}
