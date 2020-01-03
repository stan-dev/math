#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(ProbStdNormal, ccdf_log_matches_lccdf) {
  double y = 0.8;

  EXPECT_FLOAT_EQ((stan::math::std_normal_lccdf(y)),
                  (stan::math::std_normal_ccdf_log(y)));
  EXPECT_FLOAT_EQ((stan::math::std_normal_lccdf<double>(y)),
                  (stan::math::std_normal_ccdf_log<double>(y)));
}
