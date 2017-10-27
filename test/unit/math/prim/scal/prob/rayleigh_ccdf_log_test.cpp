#include <gtest/gtest.h>
#include <stan/math/prim/scal.hpp>

TEST(ProbRayleigh, ccdf_log_matches_lccdf) {
  double y = 0.8;
  double sigma = 2.3;

  EXPECT_FLOAT_EQ((stan::math::rayleigh_lccdf(y, sigma)),
                  (stan::math::rayleigh_ccdf_log(y, sigma)));
  EXPECT_FLOAT_EQ((stan::math::rayleigh_lccdf<double, double>(y, sigma)),
                  (stan::math::rayleigh_ccdf_log<double, double>(y, sigma)));
}
