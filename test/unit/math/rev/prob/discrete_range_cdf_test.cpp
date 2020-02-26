#include <stan/math/rev.hpp>
#include <gtest/gtest.h>

TEST(RevProbDiscreteRange, cdf_boundaries) {
  using stan::math::discrete_range_cdf;
  stan::math::var lower(1);
  stan::math::var upper(5);

  EXPECT_FLOAT_EQ(1.0 / 5, discrete_range_cdf(lower, lower, upper));
  EXPECT_FLOAT_EQ(1, discrete_range_cdf(upper, lower, upper));
}

TEST(RevProbDiscreteRange, out_of_support) {
  using stan::math::discrete_range_cdf;
  stan::math::var lower(1);
  stan::math::var upper(5);

  EXPECT_FLOAT_EQ(0, discrete_range_cdf(lower - 1, lower, upper));
  EXPECT_FLOAT_EQ(1, discrete_range_cdf(upper + 1, lower, upper));
}
