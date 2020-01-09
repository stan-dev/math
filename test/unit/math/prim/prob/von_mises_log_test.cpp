#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(ProbVonMises, log_matches_lpdf) {
  double y = -0.8;
  double mu = 0.4;
  double kappa = 2.3;

  EXPECT_FLOAT_EQ((stan::math::von_mises_lpdf(y, mu, kappa)),
                  (stan::math::von_mises_log(y, mu, kappa)));
  EXPECT_FLOAT_EQ((stan::math::von_mises_lpdf<true>(y, mu, kappa)),
                  (stan::math::von_mises_log<true>(y, mu, kappa)));
  EXPECT_FLOAT_EQ((stan::math::von_mises_lpdf<false>(y, mu, kappa)),
                  (stan::math::von_mises_log<false>(y, mu, kappa)));
  EXPECT_FLOAT_EQ(
      (stan::math::von_mises_lpdf<true, double, double, double>(y, mu, kappa)),
      (stan::math::von_mises_log<true, double, double, double>(y, mu, kappa)));
  EXPECT_FLOAT_EQ(
      (stan::math::von_mises_lpdf<false, double, double, double>(y, mu, kappa)),
      (stan::math::von_mises_log<false, double, double, double>(y, mu, kappa)));
  EXPECT_FLOAT_EQ(
      (stan::math::von_mises_lpdf<double, double, double>(y, mu, kappa)),
      (stan::math::von_mises_log<double, double, double>(y, mu, kappa)));
}

TEST(ProbVonMises, kappa_zero) {
  double y = -0.8;
  double mu = 0.4;
  double kappa = 0;

  EXPECT_FLOAT_EQ(-1.8378770664093453390819,
                  stan::math::von_mises_lpdf(y, mu, kappa));
}
