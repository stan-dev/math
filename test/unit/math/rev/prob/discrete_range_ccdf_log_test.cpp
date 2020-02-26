#include <stan/math/rev.hpp>
#include <gtest/gtest.h>

TEST(RevProbDiscreteRange, lccdf_boundaries) {
  using stan::math::discrete_range_lccdf;
  stan::math::var lower(1);
  stan::math::var upper(5);

  EXPECT_FLOAT_EQ(std::log(4.0 / 5), discrete_range_lccdf(lower, lower, upper));
  EXPECT_FLOAT_EQ(stan::math::LOG_ZERO,
                  discrete_range_lccdf(upper, lower, upper));
}

TEST(RevProbDiscreteRange, lccdf_out_of_support) {
  using stan::math::discrete_range_lccdf;
  stan::math::var lower(1);
  stan::math::var upper(5);

  EXPECT_FLOAT_EQ(0, discrete_range_lccdf(lower - 1, lower, upper));
  EXPECT_FLOAT_EQ(stan::math::LOG_ZERO,
                  discrete_range_lccdf(upper + 1, lower, upper));
}
