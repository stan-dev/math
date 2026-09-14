#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>

namespace {

void expect_stable_value(double actual, double expected) {
  EXPECT_NEAR(actual, expected, 1e-10 * std::max(1.0, std::abs(expected)));
}

}  // namespace

TEST_F(AgradRev, exp_mod_normal_lpdf_infinite_observations) {
  EXPECT_EQ(stan::math::NEGATIVE_INFTY,
            stan::math::exp_mod_normal_lpdf(
                stan::math::NEGATIVE_INFTY, 0.0, 1.0, 1.0));
  EXPECT_EQ(stan::math::NEGATIVE_INFTY,
            stan::math::exp_mod_normal_lpdf(stan::math::INFTY, 0.0, 1.0,
                                            1.0));
}

TEST_F(AgradRev, exp_mod_normal_lpdf_left_tail) {
  stan::math::var y = -37.537319574666846;
  stan::math::var mu = 0.0;
  stan::math::var sigma = 1.0;
  stan::math::var lambda = 1.0;

  stan::math::var logp = stan::math::exp_mod_normal_lpdf(y, mu, sigma, lambda);
  logp.grad();

  expect_stable_value(logp.val(), -709.09641828418876);
  expect_stable_value(y.adj(), 37.563233619402858);
  expect_stable_value(mu.adj(), -37.563233619402858);
  expect_stable_value(sigma.adj(), 1409.9971905846587);
  expect_stable_value(lambda.adj(), 0.97408595526398756);
}

TEST_F(AgradRev, exp_mod_normal_cdf_cancellation) {
  stan::math::var y = 0.093719047083121884;
  stan::math::var mu = 0.0;
  stan::math::var sigma = 1.0;
  stan::math::var lambda = 23.515615900485439;

  stan::math::var cdf = stan::math::exp_mod_normal_cdf(y, mu, sigma, lambda);
  cdf.grad();

  expect_stable_value(cdf.val(), 0.5204063368178875);
  expect_stable_value(y.adj(), 0.39806043386589307);
  expect_stable_value(mu.adj(), -0.39806043386589307);
  expect_stable_value(sigma.adj(), -0.020372021355659669);
  expect_stable_value(lambda.adj(), 0.000720109703246017);
}

TEST_F(AgradRev, exp_mod_normal_lcdf_cancellation) {
  stan::math::var y = 0.093719047083121884;
  stan::math::var mu = 0.0;
  stan::math::var sigma = 1.0;
  stan::math::var lambda = 23.515615900485439;

  stan::math::var lcdf = stan::math::exp_mod_normal_lcdf(y, mu, sigma, lambda);
  lcdf.grad();

  expect_stable_value(lcdf.val(), -0.65314535559646436);
  expect_stable_value(y.adj(), 0.76490312608393829);
  expect_stable_value(mu.adj(), -0.76490312608393829);
  expect_stable_value(sigma.adj(), -0.039146374504637732);
  expect_stable_value(lambda.adj(), 0.0013837450705332481);
}

TEST_F(AgradRev, exp_mod_normal_lcdf_small_lambda_lower_tail) {
  stan::math::var y = -10.0;
  stan::math::var mu = 0.0;
  stan::math::var sigma = 1.0;
  stan::math::var lambda = 1e-16;

  stan::math::var lcdf = stan::math::exp_mod_normal_lcdf(y, mu, sigma, lambda);
  lcdf.grad();

  expect_stable_value(lcdf.val(), -92.394483524027308);
  expect_stable_value(y.adj(), 10.194383033414763);
  expect_stable_value(mu.adj(), -10.194383033414763);
  expect_stable_value(sigma.adj(), 102.94383033414763);
  expect_stable_value(lambda.adj(), 1e16);
}

TEST_F(AgradRev, exp_mod_normal_lccdf_cancellation) {
  stan::math::var y = 6.2714221555374134;
  stan::math::var mu = 0.0;
  stan::math::var sigma = 1.0;
  stan::math::var lambda = 17.815654979555671;

  stan::math::var lccdf
      = stan::math::exp_mod_normal_lccdf(y, mu, sigma, lambda);
  lccdf.grad();

  expect_stable_value(lccdf.val(), -22.004519950962418);
  expect_stable_value(y.adj(), -6.339102988738702);
  expect_stable_value(mu.adj(), 6.339102988738702);
  expect_stable_value(sigma.adj(), 39.214024191656705);
  expect_stable_value(lambda.adj(), -0.030375910331314884);
}
