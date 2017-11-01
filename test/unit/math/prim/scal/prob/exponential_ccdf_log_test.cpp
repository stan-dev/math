#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(ProbExponential, ccdf_log_matches_lccdf) {
  double y = 0.8;
  double beta = 1.2;

  EXPECT_FLOAT_EQ((stan::math::exponential_lccdf(y, beta)),
                  (stan::math::exponential_ccdf_log(y, beta)));
  EXPECT_FLOAT_EQ((stan::math::exponential_lccdf<double, double>(y, beta)),
                  (stan::math::exponential_ccdf_log<double, double>(y, beta)));
}
