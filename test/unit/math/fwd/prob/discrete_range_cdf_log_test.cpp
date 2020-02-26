#include <stan/math/fwd.hpp>
#include <gtest/gtest.h>

TEST(FwdProbDiscreteRange, lcdf_boundaries) {
  using stan::math::discrete_range_lcdf;
  stan::math::fvar<double> lower(1);
  stan::math::fvar<double> upper(5);

  EXPECT_FLOAT_EQ(std::log(1.0 / 5), discrete_range_lcdf(lower, lower, upper));
  EXPECT_FLOAT_EQ(0, discrete_range_lcdf(upper, lower, upper));
}

TEST(FwdProbDiscreteRange, lcdf_out_of_support) {
  using stan::math::discrete_range_lcdf;
  stan::math::fvar<double> lower(1);
  stan::math::fvar<double> upper(5);

  EXPECT_FLOAT_EQ(stan::math::LOG_ZERO,
                  discrete_range_lcdf(lower - 1, lower, upper));
  EXPECT_FLOAT_EQ(0, discrete_range_lcdf(upper + 1, lower, upper));
}
