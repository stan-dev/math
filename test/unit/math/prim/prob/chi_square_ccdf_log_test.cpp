#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(ProbChiSquare, ccdf_log_matches_lccdf) {
  double y = 0.3;
  double nu = 15;

  EXPECT_FLOAT_EQ((stan::math::chi_square_lccdf(y, nu)),
                  (stan::math::chi_square_ccdf_log(y, nu)));
  EXPECT_FLOAT_EQ((stan::math::chi_square_lccdf<double, double>(y, nu)),
                  (stan::math::chi_square_ccdf_log<double, double>(y, nu)));
}
