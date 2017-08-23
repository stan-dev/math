#include <stan/math/prim/scal.hpp>
#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

using Eigen::Dynamic;
using Eigen::Matrix;

TEST(ProbDistributionsBernoulliLogitGLM, glm_matches_bernoulli_logit) {
  Matrix<int,Dynamic,1> n(3,1);
  n << 1, 0, 1;
  Matrix<double,Dynamic,1> x(3,1);
  x << -12, 46, -42;
  Matrix<double,1,1> beta(1,1);
  x << 0.3;
  Matrix<double,Dynamic,1> alpha(3,1);
  alpha << 10, 23, 13;
  Matrix<double,Dynamic,1> theta(3,1);
  theta = x * beta + alpha;

  EXPECT_FLOAT_EQ((stan::math::bernoulli_logit_lpmf(n, theta)),
                  (stan::math::bernoulli_logit_glm_lpmf(n, x, beta, alpha)));
  EXPECT_FLOAT_EQ((stan::math::bernoulli_logit_lpmf<true>(n, theta)),
                  (stan::math::bernoulli_logit_glm_lpmf<true>(n, x, beta, alpha)));
  EXPECT_FLOAT_EQ((stan::math::bernoulli_logit_lpmf<false>(n, theta)),
                  (stan::math::bernoulli_logit_glm_lpmf<false>(n, x, beta, alpha)));
  EXPECT_FLOAT_EQ((stan::math::bernoulli_logit_lpmf<true, double>(n, theta)),
                  (stan::math::bernoulli_logit_glm_lpmf<true, double, double, double>(n, x, beta, alpha)));
  EXPECT_FLOAT_EQ((stan::math::bernoulli_logit_lpmf<false, double>(n, theta)),
                  (stan::math::bernoulli_logit_glm_lpmf<false,  double, double, double>(n, x, beta, alpha)));
  EXPECT_FLOAT_EQ((stan::math::bernoulli_logit_lpmf<double>(n, theta)),
                  (stan::math::bernoulli_logit_glm_lpmf< double, double, double>(n, x, beta, alpha)));
}
