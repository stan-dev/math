#include <stan/math/fwd.hpp>
#include <gtest/gtest.h>

TEST(FwdProbDiscreteRange, cdf_boundaries) {
  using stan::math::discrete_range_cdf;
  stan::math::fvar<double> lower(1);
  stan::math::fvar<double> upper(5);

  EXPECT_FLOAT_EQ(1.0 / 5, discrete_range_cdf(lower, lower, upper));
  EXPECT_FLOAT_EQ(1, discrete_range_cdf(upper, lower, upper));
}

TEST(FwdProbDiscreteRange, out_of_support) {
  using stan::math::discrete_range_cdf;
  stan::math::fvar<double> lower(1);
  stan::math::fvar<double> upper(5);

  EXPECT_FLOAT_EQ(0, discrete_range_cdf(lower - 1, lower, upper));
  EXPECT_FLOAT_EQ(1, discrete_range_cdf(upper + 1, lower, upper));
}
