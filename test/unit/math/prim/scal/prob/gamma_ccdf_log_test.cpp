#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(ProbGamma, ccdf_log_matches_lccdf) {
  double y = 0.8;
  double alpha = 1.1;
  double beta = 2.3;

  EXPECT_FLOAT_EQ((stan::math::gamma_lccdf(y, alpha, beta)),
                  (stan::math::gamma_ccdf_log(y, alpha, beta)));
  EXPECT_FLOAT_EQ(
      (stan::math::gamma_lccdf<double, double, double>(y, alpha, beta)),
      (stan::math::gamma_ccdf_log<double, double, double>(y, alpha, beta)));
}
