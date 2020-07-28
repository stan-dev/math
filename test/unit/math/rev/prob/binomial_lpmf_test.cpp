#include <stan/math/rev.hpp>
#include <gtest/gtest.h>

TEST(ProbBinomial, N_equals_0) {
  using stan::math::var;
  for (double theta_val : {0.0, 0.5, 1.0}) {
    var theta = theta_val;

    var logp = stan::math::binomial_lpmf(0, 0, theta);
    logp.grad();
    EXPECT_FLOAT_EQ(logp.val(), 0.0);
    EXPECT_FLOAT_EQ(theta.adj(), 0.0);
  }

  stan::math::recover_memory();
}

TEST(ProbBinomial, N_equals_1) {
  using stan::math::var;
  for (double theta_val : {0.0, 0.5, 1.0}) {
    for (int n : {0, 1}) {
      var theta = theta_val;

      var logp = stan::math::binomial_lpmf(n, 1, theta);
      logp.grad();
      if (n == 0 && theta_val == 0.0) {
        EXPECT_FLOAT_EQ(logp.val(), 0.0);
        EXPECT_FLOAT_EQ(theta.adj(), -1.0);
      } else if (n == 1 && theta_val == 0.0) {
        EXPECT_TRUE(stan::math::is_inf(logp.val()));
      } else if (n == 0 && theta_val == 0.5) {
        EXPECT_FLOAT_EQ(logp.val(), std::log(0.5));
        EXPECT_FLOAT_EQ(theta.adj(), -2.0);
      } else if (n == 1 && theta_val == 0.5) {
        EXPECT_FLOAT_EQ(logp.val(), std::log(0.5));
        EXPECT_FLOAT_EQ(theta.adj(), 2.0);
      } else if (n == 0 && theta_val == 1.0) {
        EXPECT_TRUE(stan::math::is_inf(logp.val()));
      } else if (n == 1 && theta_val == 1.0) {
        EXPECT_FLOAT_EQ(logp.val(), 0.0);
        EXPECT_FLOAT_EQ(theta.adj(), 1.0);
      }

      stan::math::recover_memory();
    }
  }
}

TEST(ProbBinomial, N_equals_2) {
  using stan::math::var;
  for (double theta_val : {0.0, 0.5, 1.0}) {
    for (int n : {0, 1, 2}) {
      var theta = theta_val;

      var logp = stan::math::binomial_lpmf(n, 2, theta);
      logp.grad();
      if (n == 0 && theta_val == 0.0) {
        EXPECT_FLOAT_EQ(logp.val(), 0.0);
        EXPECT_FLOAT_EQ(theta.adj(), -2.0);
      } else if (n == 1 && theta_val == 0.0) {
        EXPECT_TRUE(stan::math::is_inf(logp.val()));
      } else if (n == 2 && theta_val == 0.0) {
        EXPECT_TRUE(stan::math::is_inf(logp.val()));
      } else if (n == 0 && theta_val == 0.5) {
        EXPECT_FLOAT_EQ(logp.val(), std::log(0.25));
        EXPECT_FLOAT_EQ(theta.adj(), -4.0);
      } else if (n == 1 && theta_val == 0.5) {
        EXPECT_FLOAT_EQ(logp.val(), std::log(0.5));
        EXPECT_FLOAT_EQ(theta.adj(), 0.0);
      } else if (n == 2 && theta_val == 0.5) {
        EXPECT_FLOAT_EQ(logp.val(), std::log(0.25));
        EXPECT_FLOAT_EQ(theta.adj(), 4.0);
      } else if (n == 0 && theta_val == 1.0) {
        EXPECT_TRUE(stan::math::is_inf(logp.val()));
      } else if (n == 1 && theta_val == 1.0) {
        EXPECT_TRUE(stan::math::is_inf(logp.val()));
      } else if (n == 2 && theta_val == 1.0) {
        EXPECT_FLOAT_EQ(logp.val(), 0.0);
        EXPECT_FLOAT_EQ(theta.adj(), 2.0);
      }

      stan::math::recover_memory();
    }
  }
}

TEST(ProbBinomial, n_equals_N) {
  using stan::math::var;
  var theta = 1.0;

  var logp = stan::math::binomial_lpmf(2, 2, theta);
  logp.grad();
  EXPECT_FLOAT_EQ(logp.val(), 0.0);
  EXPECT_FLOAT_EQ(theta.adj(), 2.0);

  stan::math::recover_memory();
}

TEST(ProbBinomial, n_equals_N_vec) {
  using stan::math::var;
  var theta = 1.0;
  std::vector<int> n = {2, 3};
  std::vector<int> N = {2, 3};

  var logp = stan::math::binomial_lpmf(n, N, theta);
  logp.grad();
  EXPECT_FLOAT_EQ(logp.val(), 0.0);
  EXPECT_FLOAT_EQ(theta.adj(), 5.0);

  stan::math::recover_memory();
}

TEST(ProbBinomial, n_equals_zero) {
  using stan::math::var;
  var theta = 0.0;

  var logp = stan::math::binomial_lpmf(0, 2, theta);
  logp.grad();
  EXPECT_FLOAT_EQ(logp.val(), 0.0);
  EXPECT_FLOAT_EQ(theta.adj(), -2.0);

  stan::math::recover_memory();
}

TEST(ProbBinomial, n_equals_0_vec) {
  using stan::math::var;
  var theta = 0.0;
  std::vector<int> n = {0, 0};
  std::vector<int> N = {2, 3};

  var logp = stan::math::binomial_lpmf(n, N, theta);
  logp.grad();
  EXPECT_FLOAT_EQ(logp.val(), 0.0);
  EXPECT_FLOAT_EQ(theta.adj(), -5.0);

  stan::math::recover_memory();
}

TEST(ProbBinomial, N_equals_0_vec) {
  using stan::math::var;
  var theta = 0.0;
  std::vector<int> n = {0, 0};
  std::vector<int> N = {0, 0};

  var logp = stan::math::binomial_lpmf(n, N, theta);
  logp.grad();
  EXPECT_FLOAT_EQ(logp.val(), 0.0);
  EXPECT_FLOAT_EQ(theta.adj(), 0.0);

  stan::math::recover_memory();
}
