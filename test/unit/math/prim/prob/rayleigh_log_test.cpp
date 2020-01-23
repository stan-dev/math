#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(ProbRayleigh, log_matches_lpdf) {
  double y = 0.8;
  double sigma = 2.3;

  EXPECT_FLOAT_EQ((stan::math::rayleigh_lpdf(y, sigma)),
                  (stan::math::rayleigh_log(y, sigma)));
  EXPECT_FLOAT_EQ((stan::math::rayleigh_lpdf<true>(y, sigma)),
                  (stan::math::rayleigh_log<true>(y, sigma)));
  EXPECT_FLOAT_EQ((stan::math::rayleigh_lpdf<false>(y, sigma)),
                  (stan::math::rayleigh_log<false>(y, sigma)));
  EXPECT_FLOAT_EQ((stan::math::rayleigh_lpdf<true, double, double>(y, sigma)),
                  (stan::math::rayleigh_log<true, double, double>(y, sigma)));
  EXPECT_FLOAT_EQ((stan::math::rayleigh_lpdf<false, double, double>(y, sigma)),
                  (stan::math::rayleigh_log<false, double, double>(y, sigma)));
  EXPECT_FLOAT_EQ((stan::math::rayleigh_lpdf<double, double>(y, sigma)),
                  (stan::math::rayleigh_log<double, double>(y, sigma)));
}
