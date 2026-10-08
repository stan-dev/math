#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>
#include <cmath>

// The pdf / cdf quotient underflows to 0 / 0 in the lower tail. The partials
// are now built from inv_logit(-z) and only pick up the cdf itself as a
// factor, so they stay finite all the way down to the denormal range.
TEST_F(AgradRev, logistic_cdf_lower_tail) {
  stan::math::var y = -745.0;
  stan::math::var mu = 0.0;
  stan::math::var sigma = 1.0;

  stan::math::var cdf = stan::math::logistic_cdf(y, mu, sigma);
  cdf.grad();

  // dF/dy = F * (1 - F) / sigma, and 1 - F == 1 to machine precision here
  EXPECT_GT(cdf.val(), 0.0);
  EXPECT_DOUBLE_EQ(cdf.val(), y.adj());
  EXPECT_DOUBLE_EQ(-cdf.val(), mu.adj());
  EXPECT_DOUBLE_EQ(745.0 * cdf.val(), sigma.adj());
}

// One ulp further out the cdf itself underflows to zero; the gradient must
// follow it to zero rather than become NaN.
TEST_F(AgradRev, logistic_cdf_lower_tail_underflow) {
  stan::math::var y = -746.0;
  stan::math::var mu = 0.0;
  stan::math::var sigma = 1.0;

  stan::math::var cdf = stan::math::logistic_cdf(y, mu, sigma);
  cdf.grad();

  EXPECT_DOUBLE_EQ(0.0, cdf.val());
  EXPECT_DOUBLE_EQ(0.0, y.adj());
  EXPECT_DOUBLE_EQ(0.0, mu.adj());
  EXPECT_DOUBLE_EQ(0.0, sigma.adj());
}

// inv_logit(-z) underflows to zero above z = 745; with the cdf itself equal
// to one there, the gradient is the only thing left to get right.
TEST_F(AgradRev, logistic_cdf_underflow_rescued_by_small_sigma) {
  stan::math::var y = 8e-298;
  stan::math::var mu = 0.0;
  stan::math::var sigma = 1e-300;

  stan::math::var cdf = stan::math::logistic_cdf(y, mu, sigma);
  cdf.grad();

  const double deriv = 3.6678745841780173e-48;
  EXPECT_EQ(1.0, cdf.val());
  EXPECT_NEAR(deriv, y.adj(), 1e-10 * deriv);
  EXPECT_NEAR(-deriv, mu.adj(), 1e-10 * deriv);
  EXPECT_NEAR(-2.9342996673424137e-45, sigma.adj(),
              1e-10 * 2.9342996673424137e-45);
}

TEST_F(AgradRev, logistic_cdf_location_scale) {
  const double mu_value = 5.0;
  const double sigma_value = 2.0;
  const double scaled_diff = -1.5;

  stan::math::var y = mu_value + sigma_value * scaled_diff;
  stan::math::var mu = mu_value;
  stan::math::var sigma = sigma_value;

  stan::math::var cdf = stan::math::logistic_cdf(y, mu, sigma);
  cdf.grad();

  const double cdf_value = stan::math::inv_logit(scaled_diff);
  const double deriv = stan::math::inv_logit(-scaled_diff) / sigma_value;
  EXPECT_DOUBLE_EQ(cdf_value, cdf.val());
  EXPECT_DOUBLE_EQ(deriv * cdf_value, y.adj());
  EXPECT_DOUBLE_EQ(-deriv * cdf_value, mu.adj());
  EXPECT_DOUBLE_EQ(-scaled_diff * deriv * cdf_value, sigma.adj());
}
