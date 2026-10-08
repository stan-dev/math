#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>
#include <cmath>

TEST_F(AgradRev, logistic_lpdf_location_gradient_large_values) {
  stan::math::var y = 711.0;
  stan::math::var mu = 710.0;
  stan::math::var sigma = 1.0;

  stan::math::var logp = stan::math::logistic_lpdf(y, mu, sigma);
  logp.grad();

  // z = 1: d/dmu = tanh(z / 2) / sigma, d/dy = -d/dmu,
  // d/dsigma = (z * tanh(z / 2) - 1) / sigma
  const double tanh_half = std::tanh(0.5);
  EXPECT_DOUBLE_EQ(-1.0 - 2.0 * stan::math::log1p_exp(-1.0), logp.val());
  EXPECT_DOUBLE_EQ(tanh_half, mu.adj());
  EXPECT_DOUBLE_EQ(-tanh_half, y.adj());
  EXPECT_DOUBLE_EQ(tanh_half - 1.0, sigma.adj());
}

// The 2 / (1 + exp(z)) - 1 form of the y partial returns 0 for z below eps,
// while the location partial is computed from tanh(z / 2). The two must agree
// to the last bit: they are the same quantity with opposite sign.
TEST_F(AgradRev, logistic_lpdf_gradients_near_location) {
  for (double scaled_diff : {1e-8, 1e-12, 1e-14, 1e-16}) {
    stan::math::var y = scaled_diff;
    stan::math::var mu = 0.0;
    stan::math::var sigma = 1.0;

    stan::math::var logp = stan::math::logistic_lpdf(y, mu, sigma);
    logp.grad();

    EXPECT_DOUBLE_EQ(std::tanh(0.5 * scaled_diff), mu.adj());
    EXPECT_DOUBLE_EQ(-mu.adj(), y.adj());
    stan::math::recover_memory();
  }
}

// The location partial is assigned through the container path of the edge, so
// exercise it with a vector argument and a non-unit scale.
TEST_F(AgradRev, logistic_lpdf_gradients_vectorized) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> y(2);
  y << 711.0, 712.0;
  stan::math::var mu = 710.0;
  stan::math::var sigma = 2.0;

  stan::math::var logp = stan::math::logistic_lpdf(y, mu, sigma);
  logp.grad();

  const double d0 = std::tanh(0.25) / 2.0;  // z = 0.5
  const double d1 = std::tanh(0.5) / 2.0;   // z = 1.0
  EXPECT_TRUE(std::isfinite(logp.val()));
  EXPECT_DOUBLE_EQ(d0 + d1, mu.adj());
  EXPECT_DOUBLE_EQ(-d0, y(0).adj());
  EXPECT_DOUBLE_EQ(-d1, y(1).adj());
  EXPECT_DOUBLE_EQ(
      (0.5 * std::tanh(0.25) - 1.0) / 2.0 + (1.0 * std::tanh(0.5) - 1.0) / 2.0,
      sigma.adj());
}
