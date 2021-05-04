#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(ProbDiscreteRange, cdf_boundaries) {
  using stan::math::discrete_range_cdf;
  int lower = 1;
  int upper = 5;

  EXPECT_FLOAT_EQ(1.0 / 5, discrete_range_cdf(lower, lower, upper));
  EXPECT_FLOAT_EQ(5.0 / 5, discrete_range_cdf(upper, lower, upper));
}

TEST(ProbDiscreteRange, cdf_out_of_support) {
  using stan::math::discrete_range_cdf;
  int lower = 1;
  int upper = 5;

  EXPECT_FLOAT_EQ(0.0, discrete_range_cdf(lower - 1, lower, upper));
  EXPECT_FLOAT_EQ(1.0, discrete_range_cdf(upper + 1, lower, upper));
}
