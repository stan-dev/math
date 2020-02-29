#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(ProbDiscreteRange, cdf_log_matches_lccdf) {
  using stan::math::discrete_range_ccdf_log;
  using stan::math::discrete_range_lccdf;

  for (int lower = 0; lower < 5; ++lower) {
    for (int upper = lower; upper < 5; ++upper) {
      for (int y = lower; y <= upper; ++y) {
        EXPECT_FLOAT_EQ((discrete_range_lccdf(y, lower, upper)),
                        (discrete_range_ccdf_log(y, lower, upper)));
        EXPECT_FLOAT_EQ((discrete_range_lccdf<double>(y, lower, upper)),
                        (discrete_range_ccdf_log<double>(y, lower, upper)));
      }
    }
  }
}

TEST(ProbDiscreteRange, lccdf_boundaries) {
  using stan::math::discrete_range_lccdf;
  int lower = 1;
  int upper = 5;

  EXPECT_FLOAT_EQ(std::log(4.0 / 5), discrete_range_lccdf(lower, lower, upper));
  EXPECT_FLOAT_EQ(stan::math::LOG_ZERO,
                  discrete_range_lccdf(upper, lower, upper));
}

TEST(ProbDiscreteRange, lccdf_out_of_support) {
  using stan::math::discrete_range_lccdf;
  int lower = 1;
  int upper = 5;

  EXPECT_FLOAT_EQ(0.0, discrete_range_lccdf(lower - 1, lower, upper));
  EXPECT_FLOAT_EQ(stan::math::LOG_ZERO,
                  discrete_range_lccdf(upper + 1, lower, upper));
}
