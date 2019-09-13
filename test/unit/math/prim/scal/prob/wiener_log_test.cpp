#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(ProbWiener, log_matches_lpdf) {
  double y = 2.4;
  double alpha = 2;
  double tau = 1.2;
  double beta = 0.3;
  double delta = -5;

  EXPECT_FLOAT_EQ((stan::math::wiener_lpdf(y, alpha, tau, beta, delta)),
                  (stan::math::wiener_log(y, alpha, tau, beta, delta)));
  EXPECT_FLOAT_EQ((stan::math::wiener_lpdf<true>(y, alpha, tau, beta, delta)),
                  (stan::math::wiener_log<true>(y, alpha, tau, beta, delta)));
  EXPECT_FLOAT_EQ((stan::math::wiener_lpdf<false>(y, alpha, tau, beta, delta)),
                  (stan::math::wiener_log<false>(y, alpha, tau, beta, delta)));
  EXPECT_FLOAT_EQ(
      (stan::math::wiener_lpdf<true>(
          y, alpha, tau, beta, delta)),
      (stan::math::wiener_log<true>(
          y, alpha, tau, beta, delta)));
  EXPECT_FLOAT_EQ(
      (stan::math::wiener_lpdf<false>(
          y, alpha, tau, beta, delta)),
      (stan::math::wiener_log<false>(
          y, alpha, tau, beta, delta)));
  EXPECT_FLOAT_EQ(
      (stan::math::wiener_lpdf(
          y, alpha, tau, beta, delta)),
      (stan::math::wiener_log(
          y, alpha, tau, beta, delta)));
}
