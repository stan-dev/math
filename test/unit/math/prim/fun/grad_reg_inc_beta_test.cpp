#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(grad_reg_inc_beta, 1) {
  double alpha = 1.0;
  double beta = 1.0;
  double y = 1.0;
  double digamma_alpha = stan::math::digamma(alpha);
  double digamma_beta = stan::math::digamma(beta);
  double digamma_sum = stan::math::digamma(alpha + beta);
  double betafunc = std::exp(stan::math::lbeta(alpha, beta));

  double g1 = 0;
  double g2 = 0;
  stan::math::grad_reg_inc_beta(g1, g2, alpha, beta, y, digamma_alpha,
                                digamma_beta, digamma_sum, betafunc);
  EXPECT_FLOAT_EQ(0, g1);
  EXPECT_FLOAT_EQ(-std::numeric_limits<double>::infinity(), g2);
}

TEST(grad_reg_inc_beta, 2) {
  double alpha = 1.0;
  double beta = 1.0;
  double y = 0.4;
  double digamma_alpha = stan::math::digamma(alpha);
  double digamma_beta = stan::math::digamma(beta);
  double digamma_sum = stan::math::digamma(alpha + beta);
  double betafunc = std::exp(stan::math::lbeta(alpha, beta));

  double g1 = 0;
  double g2 = 0;
  stan::math::grad_reg_inc_beta(g1, g2, alpha, beta, y, digamma_alpha,
                                digamma_beta, digamma_sum, betafunc);
  EXPECT_NEAR(-0.36651629, g1, 1e-6);
  EXPECT_NEAR(0.30649537, g2, 1e-6);
}
