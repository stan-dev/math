#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

TEST(MathFunctions, log_mix_exceptions) {
  using stan::math::log_mix;
  using stan::math::row_vector_d;
  using stan::math::vector_d;
  EXPECT_THROW(log_mix(-1, 10, 20), std::domain_error);
  EXPECT_THROW(log_mix(std::numeric_limits<double>::quiet_NaN(), 10, 20),
               std::domain_error);
  EXPECT_THROW(log_mix(0.5, std::numeric_limits<double>::quiet_NaN(), 10),
               std::domain_error);
  EXPECT_THROW(log_mix(0.5, 10, std::numeric_limits<double>::quiet_NaN()),
               std::domain_error);
}
void test_log_mix(double theta, double lambda1, double lambda2) {
  using stan::math::log_mix;
  using stan::math::row_vector_d;
  using stan::math::vector_d;
  using std::exp;
  using std::log;
  EXPECT_FLOAT_EQ(log(theta * exp(lambda1) + (1 - theta) * exp(lambda2)),
                  log_mix(theta, lambda1, lambda2));
}

TEST(MathFunctions, log_mix_values) {
  test_log_mix(0.3, 1.7, -3.9);
  test_log_mix(0.0001, 197, -3000);
  test_log_mix(0.999999, 197, -3000);
}

template <typename T_a, typename T_b>
void log_mix_val_test(T_a a, T_b b) {
  using stan::math::log_mix;
  a[0] = 0.321;
  a[1] = 0.115;
  a[2] = 0.261;
  a[3] = 0.303;

  b[0] = -5.918;
  b[1] = -7.215;
  b[2] = -9.635;
  b[3] = -8.264;

  double out = log_mix(a, b);

  EXPECT_FLOAT_EQ(out, -6.86528599744793);

  T_b b2(4), b3(4);

  b2[0] = -7.514;
  b2[1] = -2.653;
  b2[2] = -6.527;
  b2[3] = -5.618;

  b3[0] = -4.517;
  b3[1] = -6.359;
  b3[2] = -14.218;
  b3[3] = -8.628;

  std::vector<T_b> c{b, b2, b3};

  double std_out = log_mix(a, c);

  EXPECT_FLOAT_EQ(std_out, -17.0784535665594);
}

TEST(MatrixFunctions, LogMix_Combin) {
  /**
   * Test that all possible combinations of inputs return
   * the same result.
   */

  stan::math::vector_d vecd_prob(4);
  stan::math::vector_d vecd_dens(4);
  stan::math::row_vector_d row_vecd_prob(4);
  stan::math::row_vector_d row_vecd_dens(4);
  std::vector<double> std_prob(4);
  std::vector<double> std_dens(4);

  log_mix_val_test(vecd_prob, vecd_dens);
  log_mix_val_test(vecd_prob, row_vecd_dens);
  log_mix_val_test(vecd_prob, std_dens);

  log_mix_val_test(row_vecd_prob, vecd_dens);
  log_mix_val_test(row_vecd_prob, row_vecd_dens);
  log_mix_val_test(row_vecd_prob, std_dens);

  log_mix_val_test(std_prob, vecd_dens);
  log_mix_val_test(std_prob, row_vecd_dens);
  log_mix_val_test(std_prob, std_dens);
}

TEST(MatrixFunctions, LogMix_Values) {
  using stan::math::log_mix;
  /**
   * Test that the function is equivalent to scalar and
   * log_sum_exp implementations.
   */

  stan::math::vector_d prob(5, 1);
  prob << 0.1, 0.3, 0.25, 0.15, 0.2;

  stan::math::vector_d dens(5, 1);
  dens << -5.65, -7.62, -12.63, -55.62, -2.35;

  stan::math::vector_d prob_2(2, 1);
  prob_2 << 0.1, 0.9;

  stan::math::vector_d dens_2(2, 1);
  dens_2 << -5.65, -7.62;

  double log_mix_stan_1 = log_mix(prob, dens);
  double log_mix_sumexp
      = stan::math::log_sum_exp((stan::math::log(prob) + dens).eval());

  double log_mix_stan_2 = log_mix(prob_2, dens_2);
  double log_mix_stan_scal = log_mix(0.1, -5.65, -7.62);

  EXPECT_FLOAT_EQ(log_mix_stan_1, log_mix_sumexp);
  EXPECT_FLOAT_EQ(log_mix_stan_2, log_mix_stan_scal);
}

TEST(MatrixFunctions, LogMix_Throws) {
  using stan::math::log_mix;
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

  std::vector<stan::math::vector_d> std_dens{dens, dens, dens, dens};

  EXPECT_THROW(log_mix(prob_neg, dens), std::domain_error);
  EXPECT_THROW(log_mix(prob_neg, std_dens), std::domain_error);
  EXPECT_THROW(log_mix(prob_inf, dens), std::domain_error);
  EXPECT_THROW(log_mix(prob_inf, std_dens), std::domain_error);
  EXPECT_THROW(log_mix(prob_neg_inf, dens), std::domain_error);
  EXPECT_THROW(log_mix(prob_neg_inf, std_dens), std::domain_error);
  EXPECT_THROW(log_mix(prob_nan, dens), std::domain_error);
  EXPECT_THROW(log_mix(prob_nan, std_dens), std::domain_error);

  /**
   * Test invalid vector of densities
   */
  stan::math::vector_d dens_inf(5, 1);
  dens_inf << stan::math::INFTY, -7.62, -12.63, -55.62, -2.35;
  std::vector<stan::math::vector_d> std_dens_inf{dens_inf, dens_inf, dens_inf};

  stan::math::vector_d dens_neg_inf(5, 1);
  dens_neg_inf << -5.65, -7.62, stan::math::NEGATIVE_INFTY, -55.62, -2.35;
  std::vector<stan::math::vector_d> std_dens_neg_inf{dens_neg_inf, dens_neg_inf,
                                                     dens_neg_inf};

  stan::math::vector_d dens_nan(5, 1);
  dens_nan << -5.65, -7.62, -12.63, stan::math::NOT_A_NUMBER, -2.35;
  std::vector<stan::math::vector_d> std_dens_nan{dens_nan, dens_nan, dens_nan};

  stan::math::vector_d prob(5, 1);
  prob << 0.1, 0.3, 0.25, 0.15, 0.2;

  EXPECT_THROW(log_mix(prob, dens_inf), std::domain_error);
  EXPECT_THROW(log_mix(prob, std_dens_inf), std::domain_error);
  EXPECT_THROW(log_mix(prob, dens_neg_inf), std::domain_error);
  EXPECT_THROW(log_mix(prob, std_dens_neg_inf), std::domain_error);
  EXPECT_THROW(log_mix(prob, dens_nan), std::domain_error);
  EXPECT_THROW(log_mix(prob, std_dens_nan), std::domain_error);

  /**
   * Test inconsistent vector sizes
   */
  stan::math::vector_d prob_small(2, 1);
  prob_small << 0.5, 0.5;

  stan::math::vector_d dens_large(7, 1);
  dens_large << -5.65, -7.62, -12.63, -55.62, -2.35, -9.6, -0.7;
  std::vector<stan::math::vector_d> std_dens_large{dens_large, dens_large,
                                                   dens_large};

  EXPECT_THROW(log_mix(prob_small, dens), std::invalid_argument);
  EXPECT_THROW(log_mix(prob_small, std_dens), std::invalid_argument);
  EXPECT_THROW(log_mix(prob, dens_large), std::invalid_argument);
  EXPECT_THROW(log_mix(prob, std_dens_large), std::invalid_argument);
}
