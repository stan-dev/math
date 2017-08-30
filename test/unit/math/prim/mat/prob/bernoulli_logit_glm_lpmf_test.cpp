//#include <stan/math/prim/mat.hpp>
#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <stan/math/prim/scal/fun/value_of.hpp>

using stan::math::var;
using Eigen::Dynamic;
using Eigen::Matrix;

TEST(ProbDistributionsBernoulliLogitGLM, glm_matches_bernoulli_logit_doubles) {
  Matrix<int,Dynamic,1> n(3,1);
  n << 1, 0, 1;
  Matrix<double,Dynamic,Dynamic> x(3,2);
  x << -12, 46, -42,
       24, 25, 27; 
  Matrix<double,Dynamic,1> beta(2,1);
  beta << 0.3, 2;
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
  EXPECT_FLOAT_EQ((stan::math::bernoulli_logit_lpmf<true, Matrix<int,Dynamic,1> >(n, theta)),
                  (stan::math::bernoulli_logit_glm_lpmf<true, Matrix<int,Dynamic,1> >(n, x, beta, alpha)));
  EXPECT_FLOAT_EQ((stan::math::bernoulli_logit_lpmf<false, Matrix<int,Dynamic,1> >(n, theta)),
                  (stan::math::bernoulli_logit_glm_lpmf<false,  Matrix<int,Dynamic,1>>(n, x, beta, alpha)));
  EXPECT_FLOAT_EQ((stan::math::bernoulli_logit_lpmf<Matrix<int,Dynamic,1> >(n, theta)),
                  (stan::math::bernoulli_logit_glm_lpmf< Matrix<int,Dynamic,1> >(n, x, beta, alpha)));
				  
}

TEST(ProbDistributionsBernoulliLogitGLM, glm_matches_bernoulli_logit_vars) {
  Matrix<int,Dynamic,1> n(3,1);
  n << 1, 0, 1;
  Matrix<var,Dynamic,Dynamic> x(3,2);
  x << -12, 46, -42,
       24, 25, 27;
  Matrix<var,Dynamic,1> beta(2,1);
  beta << 0.3, 2;
  Matrix<var,Dynamic,1> alpha(3,1);
  alpha << 10, 23, 13;
  Matrix<var,Dynamic,1> theta(3,1);
  theta = x * beta + alpha;
  var lp = stan::math::bernoulli_logit_lpmf(n, theta);
  lp.grad();
  stan::math::recover_memory();
  
  Matrix<int,Dynamic,1> n2(3,1);
  n2 << 1, 0, 1;
  Matrix<var,Dynamic,Dynamic> x2(3,2);
  x2 << -12, 46, -42,
       24, 25, 27; 
  Matrix<var,Dynamic,1> beta2(2,1);
  beta2 << 0.3, 2;
  Matrix<var,Dynamic,1> alpha2(3,1);
  alpha2 << 10, 23, 13;
  var lp2 = stan::math::bernoulli_logit_glm_lpmf(n2, x2, beta2, alpha2);
  lp2.grad();
  
  std::cout << beta[0].adj();
  std::cout << beta2[0].adj();
  
  EXPECT_FLOAT_EQ(value_of(stan::math::bernoulli_logit_lpmf(n, theta)),
                  value_of(stan::math::bernoulli_logit_glm_lpmf(n, x, beta, alpha)));
  /*EXPECT_FLOAT_EQ((stan::math::bernoulli_logit_lpmf<true>(n, theta)),
                  (stan::math::bernoulli_logit_glm_lpmf<true>(n, x, beta, alpha)));
  EXPECT_FLOAT_EQ((stan::math::bernoulli_logit_lpmf<false>(n, theta)),
                  (stan::math::bernoulli_logit_glm_lpmf<false>(n, x, beta, alpha)));
  EXPECT_FLOAT_EQ((stan::math::bernoulli_logit_lpmf<true, Matrix<int,Dynamic,1> >(n, theta)),
                  (stan::math::bernoulli_logit_glm_lpmf<true, Matrix<int,Dynamic,1> >(n, x, beta, alpha)));
  EXPECT_FLOAT_EQ((stan::math::bernoulli_logit_lpmf<false, Matrix<int,Dynamic,1> >(n, theta)),
                  (stan::math::bernoulli_logit_glm_lpmf<false,  Matrix<int,Dynamic,1>>(n, x, beta, alpha)));
  EXPECT_FLOAT_EQ((stan::math::bernoulli_logit_lpmf<Matrix<int,Dynamic,1> >(n, theta)),
                  (stan::math::bernoulli_logit_glm_lpmf< Matrix<int,Dynamic,1> >(n, x, beta, alpha)));*/
				  
}
