#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>
#include <cmath>

TEST_F(AgradRev, logistic_lccdf_upper_tail) {
  for (double y_value : {36.5, 40.0}) {
    stan::math::var y = y_value;
    stan::math::var mu = 0.0;
    stan::math::var sigma = 1.0;

    stan::math::var log_ccdf = stan::math::logistic_lccdf(y, mu, sigma);
    log_ccdf.grad();

    const double deriv = stan::math::inv_logit(y_value);
    EXPECT_DOUBLE_EQ(-stan::math::log1p_exp(y_value), log_ccdf.val());
    EXPECT_DOUBLE_EQ(-deriv, y.adj());
    EXPECT_DOUBLE_EQ(deriv, mu.adj());
    EXPECT_DOUBLE_EQ(y_value * deriv, sigma.adj());
    stan::math::recover_memory();
  }
}

// Same tail, but with mu != 0 and sigma != 1 so that the scaled difference is
// distinguishable from y and the 1 / sigma factors in the partials are
// exercised.
TEST_F(AgradRev, logistic_lccdf_upper_tail_location_scale) {
  const double mu_value = 5.0;
  const double sigma_value = 2.0;

  for (double scaled_diff : {36.5, 40.0}) {
    stan::math::var y = mu_value + sigma_value * scaled_diff;
    stan::math::var mu = mu_value;
    stan::math::var sigma = sigma_value;

    stan::math::var log_ccdf = stan::math::logistic_lccdf(y, mu, sigma);
    log_ccdf.grad();

    const double deriv = stan::math::inv_logit(scaled_diff) / sigma_value;
    EXPECT_DOUBLE_EQ(-stan::math::log1p_exp(scaled_diff), log_ccdf.val());
    EXPECT_DOUBLE_EQ(-deriv, y.adj());
    EXPECT_DOUBLE_EQ(deriv, mu.adj());
    EXPECT_DOUBLE_EQ(scaled_diff * deriv, sigma.adj());
    stan::math::recover_memory();
  }
}

// The 1 - inv_logit(z) cancellation degrades long before it returns -Inf:
// at z = 30 it gave -30.001021 for the value and -1.00102 for the y partial.
TEST_F(AgradRev, logistic_lccdf_upper_tail_moderate) {
  stan::math::var y = 30.0;
  stan::math::var mu = 0.0;
  stan::math::var sigma = 1.0;

  stan::math::var log_ccdf = stan::math::logistic_lccdf(y, mu, sigma);
  log_ccdf.grad();

  EXPECT_NEAR(-30.000000000000092, log_ccdf.val(), 1e-12);
  EXPECT_NEAR(-0.99999999999990652, y.adj(), 1e-12);
  EXPECT_NEAR(0.99999999999990652, mu.adj(), 1e-12);
  EXPECT_NEAR(29.999999999997197, sigma.adj(), 1e-10);
}

// inv_logit(z) underflows to zero below z = -745, but dividing by a small
// enough sigma brings the quotient back into range.
TEST_F(AgradRev, logistic_lccdf_underflow_rescued_by_small_sigma) {
  stan::math::var y = -8e-298;
  stan::math::var mu = 0.0;
  stan::math::var sigma = 1e-300;

  stan::math::var log_ccdf = stan::math::logistic_lccdf(y, mu, sigma);
  log_ccdf.grad();

  const double deriv = 3.6678745841780173e-48;
  EXPECT_EQ(0.0, stan::math::inv_logit(-8e-298 / 1e-300));
  EXPECT_NEAR(-deriv, y.adj(), 1e-10 * deriv);
  EXPECT_NEAR(deriv, mu.adj(), 1e-10 * deriv);
  EXPECT_NEAR(-2.9342996673424137e-45, sigma.adj(),
              1e-10 * 2.9342996673424137e-45);
}

// An infinite element short-circuits the result; the partials of the finite
// elements that precede it must not leak into the returned gradient.
TEST_F(AgradRev, logistic_lccdf_pos_inf_zeroes_partials) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> y(2);
  y << 1.5, stan::math::INFTY;
  stan::math::var mu = 0.0;
  stan::math::var sigma = 1.0;

  stan::math::var log_ccdf = stan::math::logistic_lccdf(y, mu, sigma);
  log_ccdf.grad();

  EXPECT_TRUE(std::isinf(log_ccdf.val()));
  EXPECT_LT(log_ccdf.val(), 0.0);
  EXPECT_DOUBLE_EQ(0.0, y(0).adj());
  EXPECT_DOUBLE_EQ(0.0, y(1).adj());
  EXPECT_DOUBLE_EQ(0.0, mu.adj());
  EXPECT_DOUBLE_EQ(0.0, sigma.adj());
}
