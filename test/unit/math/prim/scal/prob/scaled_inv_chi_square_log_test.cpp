#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(ProbScaledInvChiSquare, log_matches_lpdf) {
  double y = 0.8;
  double nu = 1.1;
  double sigma = 2.3;

  EXPECT_FLOAT_EQ((stan::math::scaled_inv_chi_square_lpdf(y, nu, sigma)),
                  (stan::math::scaled_inv_chi_square_log(y, nu, sigma)));
  EXPECT_FLOAT_EQ((stan::math::scaled_inv_chi_square_lpdf<true>(y, nu, sigma)),
                  (stan::math::scaled_inv_chi_square_log<true>(y, nu, sigma)));
  EXPECT_FLOAT_EQ((stan::math::scaled_inv_chi_square_lpdf<false>(y, nu, sigma)),
                  (stan::math::scaled_inv_chi_square_log<false>(y, nu, sigma)));
  EXPECT_FLOAT_EQ(
      (stan::math::scaled_inv_chi_square_lpdf<true, double, double, double>(
          y, nu, sigma)),
      (stan::math::scaled_inv_chi_square_log<true, double, double, double>(
          y, nu, sigma)));
  EXPECT_FLOAT_EQ(
      (stan::math::scaled_inv_chi_square_lpdf<false, double, double, double>(
          y, nu, sigma)),
      (stan::math::scaled_inv_chi_square_log<false, double, double, double>(
          y, nu, sigma)));
  EXPECT_FLOAT_EQ(
      (stan::math::scaled_inv_chi_square_lpdf<double, double, double>(y, nu,
                                                                      sigma)),
      (stan::math::scaled_inv_chi_square_log<double, double, double>(y, nu,
                                                                     sigma)));
}
