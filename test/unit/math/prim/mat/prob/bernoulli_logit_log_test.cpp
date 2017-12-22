#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

using Eigen::Dynamic;
using Eigen::Matrix;

TEST(ProbBernoulliLogit, log_matches_lpmf) {
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
