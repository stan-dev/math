#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(ProbUniform, ccdf_log_matches_lccdf) {
  double y = 0.8;
  double alpha = 0.4;
  double beta = 2.3;

  EXPECT_FLOAT_EQ((stan::math::uniform_lccdf(y, alpha, beta)),
                  (stan::math::uniform_ccdf_log(y, alpha, beta)));
  EXPECT_FLOAT_EQ((stan::math::uniform_lccdf<double, double, double>(y, alpha, beta)),
                  (stan::math::uniform_ccdf_log<double, double, double>(y, alpha, beta)));
}
