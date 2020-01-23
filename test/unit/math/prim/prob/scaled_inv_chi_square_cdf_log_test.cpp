#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(ProbScaledInvChiSquare, cdf_log_matches_lcdf) {
  double y = 0.8;
  double nu = 1.1;
  double sigma = 2.3;

  EXPECT_FLOAT_EQ((stan::math::scaled_inv_chi_square_lcdf(y, nu, sigma)),
                  (stan::math::scaled_inv_chi_square_cdf_log(y, nu, sigma)));
  EXPECT_FLOAT_EQ(
      (stan::math::scaled_inv_chi_square_lcdf<double, double, double>(y, nu,
                                                                      sigma)),
      (stan::math::scaled_inv_chi_square_cdf_log<double, double, double>(
          y, nu, sigma)));
}
