#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>
#include <cmath>

TEST_F(AgradRev, logistic_lpdf_location_gradient_large_values) {
  stan::math::var mu = 710.0;

  stan::math::var logp = stan::math::logistic_lpdf(711.0, mu, 1.0);
  logp.grad();

  EXPECT_TRUE(std::isfinite(logp.val()));
  EXPECT_DOUBLE_EQ(std::tanh(0.5), mu.adj());
}
