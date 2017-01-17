#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(ProbUniform, cdf_log_matches_lcdf) {
  double y = 0.8;
  double alpha = 0.4;
  double beta = 2.3;

  EXPECT_FLOAT_EQ((stan::math::uniform_lcdf(y, alpha, beta)),
                  (stan::math::uniform_cdf_log(y, alpha, beta)));
  EXPECT_FLOAT_EQ((stan::math::uniform_lcdf<double, double, double>(y, alpha, beta)),
                  (stan::math::uniform_cdf_log<double, double, double>(y, alpha, beta)));
}
