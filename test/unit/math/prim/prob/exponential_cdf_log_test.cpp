#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(ProbExponential, cdf_log_matches_lcdf) {
  double y = 0.8;
  double beta = 1.2;

  EXPECT_FLOAT_EQ((stan::math::exponential_lcdf(y, beta)),
                  (stan::math::exponential_cdf_log(y, beta)));
  EXPECT_FLOAT_EQ((stan::math::exponential_lcdf<double, double>(y, beta)),
                  (stan::math::exponential_cdf_log<double, double>(y, beta)));
}
