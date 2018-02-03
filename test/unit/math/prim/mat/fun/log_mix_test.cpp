#include <stan/math/prim/mat/fun/log_mix.hpp>
#include <stan/math/prim/scal/fun/log_mix.hpp>
#include <stan/math/prim/mat/fun/log_sum_exp.hpp>
#include <stan/math/prim/mat/fun/log.hpp>
#include <stan/math/prim/mat/fun/typedefs.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

TEST(MatrixFunctions, LogMix_Values) {
  stan::math::vector_d prob(5, 1);
  prob << 0.1, 0.3, 0.25, 0.15, 0.2;

  std::vector<double> std_dens(5);
  std_dens[0] = -5.65;
  std_dens[1] = -7.62;
  std_dens[2] = -12.63;
  std_dens[3] = -55.62;
  std_dens[4] = -2.35;

  stan::math::vector_d dens(5, 1);
  dens << -5.65, -7.62, -12.63, -55.62, -2.35;

  stan::math::vector_d prob_2(2, 1);
  prob_2 << 0.1, 0.9;

  stan::math::vector_d dens_2(2, 1);
  dens_2 << -5.65, -7.62;

  stan::math::vector_d tmp(5, 1);
  tmp = (stan::math::log(prob) + dens);

  double log_mix_stan_1 = stan::math::log_mix(prob, dens);
  double log_mix_sumexp = stan::math::log_sum_exp(tmp);

  double log_mix_stan_2 = stan::math::log_mix(prob_2, dens_2);
  double log_mix_stan_scal = stan::math::log_mix(0.1, -5.65, -7.62);
  double log_mix_stan_3 = stan::math::log_mix(prob, std_dens);

  EXPECT_FLOAT_EQ(log_mix_stan_1, log_mix_sumexp);
  EXPECT_FLOAT_EQ(log_mix_stan_1, log_mix_stan_3);
  EXPECT_FLOAT_EQ(log_mix_stan_2, log_mix_stan_scal);
}

TEST(MatrixFunctions, LogMix_Throws) {
  /**
   * Test invalid vector of probabilities
   */
  stan::math::vector_d prob_neg(5, 1);
  prob_neg << -0.1, 0.3, 0.25, 0.15, 0.2;

  stan::math::vector_d prob_inf(5, 1);
  prob_inf << stan::math::INFTY, 0.3, 0.25, 0.15, 0.2;

  stan::math::vector_d prob_neg_inf(5, 1);
  prob_neg_inf << stan::math::NEGATIVE_INFTY, 0.3, 0.25, 0.15, 0.2;

  stan::math::vector_d prob_nan(5, 1);
  prob_nan << stan::math::NOT_A_NUMBER, 0.3, 0.25, 0.15, 0.2;

  stan::math::vector_d dens(5, 1);
  dens << -5.65, -7.62, -12.63, -55.62, -2.35;

  EXPECT_THROW(stan::math::log_mix(prob_neg, dens), std::domain_error);
  EXPECT_THROW(stan::math::log_mix(prob_inf, dens), std::domain_error);
  EXPECT_THROW(stan::math::log_mix(prob_neg_inf, dens), std::domain_error);
  EXPECT_THROW(stan::math::log_mix(prob_nan, dens), std::domain_error);

  /**
   * Test invalid vector of densities
   */
  stan::math::vector_d dens_inf(5, 1);
  dens_inf << stan::math::INFTY, -7.62, -12.63, -55.62, -2.35;

  stan::math::vector_d dens_neg_inf(5, 1);
  dens_neg_inf << -5.65, -7.62, stan::math::NEGATIVE_INFTY, -55.62, -2.35;

  stan::math::vector_d dens_nan(5, 1);
  dens_nan << -5.65, -7.62, -12.63, stan::math::NOT_A_NUMBER, -2.35;

  stan::math::vector_d prob(5, 1);
  prob << 0.1, 0.3, 0.25, 0.15, 0.2;

  EXPECT_THROW(stan::math::log_mix(prob, dens_inf), std::domain_error);
  EXPECT_THROW(stan::math::log_mix(prob, dens_neg_inf), std::domain_error);
  EXPECT_THROW(stan::math::log_mix(prob, dens_nan), std::domain_error);

  /**
   * Test inconsistent vector sizes
   */
  stan::math::vector_d prob_small(2, 1);
  prob_small << 0.5, 0.5;

  stan::math::vector_d dens_large(7, 1);
  dens_large << -5.65, -7.62, -12.63, -55.62, -2.35, -9.6, -0.7;

  EXPECT_THROW(stan::math::log_mix(prob_small, dens), std::invalid_argument);
  EXPECT_THROW(stan::math::log_mix(prob, dens_large), std::invalid_argument);
}
