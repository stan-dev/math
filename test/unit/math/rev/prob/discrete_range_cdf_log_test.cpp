#include <stan/math/rev.hpp>
#include <gtest/gtest.h>

TEST(RevProbDiscreteRange, lcdf_boundaries) {
  using stan::math::discrete_range_lcdf;
  stan::math::var lower(1);
  stan::math::var upper(5);

  EXPECT_FLOAT_EQ(std::log(1.0 / 5), discrete_range_lcdf(lower, lower, upper));
  EXPECT_FLOAT_EQ(0, discrete_range_lcdf(upper, lower, upper));
}

TEST(RevProbDiscreteRange, lcdf_out_of_support) {
  using stan::math::discrete_range_lcdf;
  stan::math::var lower(1);
  stan::math::var upper(5);

  EXPECT_FLOAT_EQ(stan::math::LOG_ZERO,
                  discrete_range_lcdf(lower - 1, lower, upper));
  EXPECT_FLOAT_EQ(0, discrete_range_lcdf(upper + 1, lower, upper));
}
