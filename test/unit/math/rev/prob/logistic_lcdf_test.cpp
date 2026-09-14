#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>

TEST_F(AgradRev, logistic_lcdf_lower_tail) {
  for (double y_value : {-745.0, -746.0}) {
    stan::math::var y = y_value;
    stan::math::var mu = 0.0;
    stan::math::var sigma = 1.0;

    stan::math::var log_cdf = stan::math::logistic_lcdf(y, mu, sigma);
    log_cdf.grad();

    const double deriv = stan::math::inv_logit(-y_value);
    EXPECT_DOUBLE_EQ(-stan::math::log1p_exp(-y_value), log_cdf.val());
    EXPECT_DOUBLE_EQ(deriv, y.adj());
    EXPECT_DOUBLE_EQ(-deriv, mu.adj());
    EXPECT_DOUBLE_EQ(-y_value * deriv, sigma.adj());
    stan::math::recover_memory();
  }
}

// Same tail, but with mu != 0 and sigma != 1 so that the scaled difference is
// distinguishable from y and the 1 / sigma factors in the partials are
// exercised.
TEST_F(AgradRev, logistic_lcdf_lower_tail_location_scale) {
  const double mu_value = 5.0;
  const double sigma_value = 2.0;

  for (double scaled_diff : {-745.0, -750.0}) {
    stan::math::var y = mu_value + sigma_value * scaled_diff;
    stan::math::var mu = mu_value;
    stan::math::var sigma = sigma_value;

    stan::math::var log_cdf = stan::math::logistic_lcdf(y, mu, sigma);
    log_cdf.grad();

    const double deriv = stan::math::inv_logit(-scaled_diff) / sigma_value;
    EXPECT_DOUBLE_EQ(-stan::math::log1p_exp(-scaled_diff), log_cdf.val());
    EXPECT_DOUBLE_EQ(deriv, y.adj());
    EXPECT_DOUBLE_EQ(-deriv, mu.adj());
    EXPECT_DOUBLE_EQ(-scaled_diff * deriv, sigma.adj());
    stan::math::recover_memory();
  }
}

// y == INFTY contributes log(1) = 0 and zero partials; it must not poison the
// gradient of the finite elements.
TEST_F(AgradRev, logistic_lcdf_pos_inf_element) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> y(2);
  y << 1.5, stan::math::INFTY;
  stan::math::var mu = 0.0;
  stan::math::var sigma = 1.0;

  stan::math::var log_cdf = stan::math::logistic_lcdf(y, mu, sigma);
  log_cdf.grad();

  const double deriv = stan::math::inv_logit(-1.5);
  EXPECT_DOUBLE_EQ(stan::math::log_inv_logit(1.5), log_cdf.val());
  EXPECT_DOUBLE_EQ(deriv, y(0).adj());
  EXPECT_DOUBLE_EQ(0.0, y(1).adj());
  EXPECT_DOUBLE_EQ(-deriv, mu.adj());
  EXPECT_DOUBLE_EQ(-1.5 * deriv, sigma.adj());
}
