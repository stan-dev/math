#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(ProbChiSquare, cdf_log_matches_lcdf) {
  double y = 0.3;
  double nu = 15;

  EXPECT_FLOAT_EQ((stan::math::chi_square_lcdf(y, nu)),
                  (stan::math::chi_square_cdf_log(y, nu)));
  EXPECT_FLOAT_EQ((stan::math::chi_square_lcdf<double, double>(y, nu)),
                  (stan::math::chi_square_cdf_log<double, double>(y, nu)));
}
