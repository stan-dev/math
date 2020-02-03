#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(ProbInvChiSquare, ccdf_log_matches_lccdf) {
  double y = 0.8;
  double nu = 2.3;

  EXPECT_FLOAT_EQ((stan::math::inv_chi_square_lccdf(y, nu)),
                  (stan::math::inv_chi_square_ccdf_log(y, nu)));
  EXPECT_FLOAT_EQ((stan::math::inv_chi_square_lccdf<double, double>(y, nu)),
                  (stan::math::inv_chi_square_ccdf_log<double, double>(y, nu)));
}
