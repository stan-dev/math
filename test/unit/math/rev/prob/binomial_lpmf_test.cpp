#include <stan/math/rev.hpp>
#include <gtest/gtest.h>

TEST(ProbBinomial, n_equals_N) {
  using stan::math::var;
  var theta = 1.0;

  var logp = stan::math::binomial_lpmf(2, 2, theta);
  logp.grad();
  EXPECT_FLOAT_EQ(logp.val(), 0.0);
  EXPECT_FLOAT_EQ(theta.adj(), 2.0);
}

TEST(ProbBinomial, n_equals_N_vec) {
  using stan::math::var;
  var theta = 1.0;
  std::vector<int> n = { 2, 3 };
  std::vector<int> N = { 2, 3 };

  var logp = stan::math::binomial_lpmf(n, N, theta);
  logp.grad();
  EXPECT_FLOAT_EQ(logp.val(), 0.0);
  EXPECT_FLOAT_EQ(theta.adj(), 5.0);
}

TEST(ProbBinomial, n_equals_zero) {
  using stan::math::var;
  var theta = 0.0;

  var logp = stan::math::binomial_lpmf(0, 2, theta);
  logp.grad();
  EXPECT_FLOAT_EQ(logp.val(), 0.0);
  EXPECT_FLOAT_EQ(theta.adj(), -2.0);
}

TEST(ProbBinomial, n_equals_0_vec) {
  using stan::math::var;
  var theta = 0.0;
  std::vector<int> n = { 0, 0 };
  std::vector<int> N = { 2, 3 };

  var logp = stan::math::binomial_lpmf(n, N, theta);
  logp.grad();
  EXPECT_FLOAT_EQ(logp.val(), 0.0);
  EXPECT_FLOAT_EQ(theta.adj(), -5.0);
}
