#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(ProbBernoulliLogitMat, log_matches_lpmf) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  Matrix<int, Dynamic, 1> n(3, 1);
  n << 0, 1, 0;
  Matrix<double, Dynamic, 1> theta(3, 1);
  theta << 1.2, 2, 0.9;

  EXPECT_FLOAT_EQ((stan::math::bernoulli_logit_lpmf(n, theta)),
                  (stan::math::bernoulli_logit_log(n, theta)));
  EXPECT_FLOAT_EQ((stan::math::bernoulli_logit_lpmf<true>(n, theta)),
                  (stan::math::bernoulli_logit_log<true>(n, theta)));
  EXPECT_FLOAT_EQ((stan::math::bernoulli_logit_lpmf<false>(n, theta)),
                  (stan::math::bernoulli_logit_log<false>(n, theta)));
  EXPECT_FLOAT_EQ(
      (stan::math::bernoulli_logit_lpmf<true, Matrix<int, Dynamic, 1>>(n,
                                                                       theta)),
      (stan::math::bernoulli_logit_log<true, Matrix<int, Dynamic, 1>>(n,
                                                                      theta)));
  EXPECT_FLOAT_EQ(
      (stan::math::bernoulli_logit_lpmf<false, Matrix<int, Dynamic, 1>>(n,
                                                                        theta)),
      (stan::math::bernoulli_logit_log<false, Matrix<int, Dynamic, 1>>(n,
                                                                       theta)));
  EXPECT_FLOAT_EQ(
      (stan::math::bernoulli_logit_lpmf<Matrix<int, Dynamic, 1>>(n, theta)),
      (stan::math::bernoulli_logit_log<Matrix<int, Dynamic, 1>>(n, theta)));
}

TEST(ProbBernoulliLogitScal, log_matches_lpmf) {
  int n = 1;
  double theta = 1.2;

  EXPECT_FLOAT_EQ((stan::math::bernoulli_logit_lpmf(n, theta)),
                  (stan::math::bernoulli_logit_log(n, theta)));
  EXPECT_FLOAT_EQ((stan::math::bernoulli_logit_lpmf<true>(n, theta)),
                  (stan::math::bernoulli_logit_log<true>(n, theta)));
  EXPECT_FLOAT_EQ((stan::math::bernoulli_logit_lpmf<false>(n, theta)),
                  (stan::math::bernoulli_logit_log<false>(n, theta)));
  EXPECT_FLOAT_EQ((stan::math::bernoulli_logit_lpmf<true, double>(n, theta)),
                  (stan::math::bernoulli_logit_log<true, double>(n, theta)));
  EXPECT_FLOAT_EQ((stan::math::bernoulli_logit_lpmf<false, double>(n, theta)),
                  (stan::math::bernoulli_logit_log<false, double>(n, theta)));
  EXPECT_FLOAT_EQ((stan::math::bernoulli_logit_lpmf<double>(n, theta)),
                  (stan::math::bernoulli_logit_log<double>(n, theta)));
}
