#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(ProbInvChiSquare, cdf_log_matches_lcdf) {
  double y = 0.8;
  double nu = 2.3;

  EXPECT_FLOAT_EQ((stan::math::inv_chi_square_lcdf(y, nu)),
                  (stan::math::inv_chi_square_cdf_log(y, nu)));
  EXPECT_FLOAT_EQ((stan::math::inv_chi_square_lcdf<double, double>(y, nu)),
                  (stan::math::inv_chi_square_cdf_log<double, double>(y, nu)));
}
