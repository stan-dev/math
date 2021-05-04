#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(ProbRayleigh, cdf_log_matches_lcdf) {
  double y = 0.8;
  double sigma = 2.3;

  EXPECT_FLOAT_EQ((stan::math::rayleigh_lcdf(y, sigma)),
                  (stan::math::rayleigh_cdf_log(y, sigma)));
  EXPECT_FLOAT_EQ((stan::math::rayleigh_lcdf<double, double>(y, sigma)),
                  (stan::math::rayleigh_cdf_log<double, double>(y, sigma)));
}
