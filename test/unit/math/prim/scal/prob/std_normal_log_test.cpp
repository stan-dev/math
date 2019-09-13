#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(ProbStdNormal, log_matches_lpdf) {
  double y = 0.8;

  EXPECT_FLOAT_EQ((stan::math::std_normal_lpdf(y)),
                  (stan::math::std_normal_log(y)));
  EXPECT_FLOAT_EQ((stan::math::std_normal_lpdf<true>(y)),
                  (stan::math::std_normal_log<true>(y)));
  EXPECT_FLOAT_EQ((stan::math::std_normal_lpdf<false>(y)),
                  (stan::math::std_normal_log<false>(y)));
  EXPECT_FLOAT_EQ((stan::math::std_normal_lpdf<true>(y)),
                  (stan::math::std_normal_log<true>(y)));
  EXPECT_FLOAT_EQ((stan::math::std_normal_lpdf<false>(y)),
                  (stan::math::std_normal_log<false>(y)));
  EXPECT_FLOAT_EQ((stan::math::std_normal_lpdf(y)),
                  (stan::math::std_normal_log(y)));
}
