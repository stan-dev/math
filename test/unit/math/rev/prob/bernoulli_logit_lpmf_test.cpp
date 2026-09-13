#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <vector>

TEST_F(AgradRev, bernoulli_logit_lpmf_upper_tail_gradient_scalar) {
  stan::math::var theta = 25.0;

  stan::math::var logp = stan::math::bernoulli_logit_lpmf(1, theta);
  logp.grad();

  EXPECT_DOUBLE_EQ(std::exp(-25.0), theta.adj());
}

TEST_F(AgradRev, bernoulli_logit_lpmf_upper_tail_gradient_vector) {
  std::vector<int> n{1, 0};
  std::vector<stan::math::var> theta{25.0, -25.0};

  stan::math::var logp = stan::math::bernoulli_logit_lpmf(n, theta);
  logp.grad();

  EXPECT_DOUBLE_EQ(std::exp(-25.0), theta[0].adj());
  EXPECT_DOUBLE_EQ(-std::exp(-25.0), theta[1].adj());
}
