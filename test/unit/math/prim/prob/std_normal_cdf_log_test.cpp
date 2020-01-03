#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(ProbStdNormal, cdf_log_matches_lcdf) {
  double y = 0.8;

  EXPECT_FLOAT_EQ((stan::math::std_normal_lcdf(y)),
                  (stan::math::std_normal_cdf_log(y)));
  EXPECT_FLOAT_EQ((stan::math::std_normal_lcdf<double>(y)),
                  (stan::math::std_normal_cdf_log<double>(y)));
}
