#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(ProbScaledInvChiSquare, ccdf_log_matches_lccdf) {
  double y = 0.8;
  double nu = 1.1;
  double sigma = 2.3;

  EXPECT_FLOAT_EQ((stan::math::scaled_inv_chi_square_lccdf(y, nu, sigma)),
                  (stan::math::scaled_inv_chi_square_ccdf_log(y, nu, sigma)));
  EXPECT_FLOAT_EQ((stan::math::scaled_inv_chi_square_lccdf<double, double, double>(y, nu, sigma)),
                  (stan::math::scaled_inv_chi_square_ccdf_log<double, double, double>(y, nu, sigma)));
}
