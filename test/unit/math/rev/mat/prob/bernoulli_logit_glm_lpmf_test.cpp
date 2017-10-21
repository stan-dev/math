#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <chrono>

using stan::math::var;
using Eigen::Dynamic;
using Eigen::Matrix;

typedef std::chrono::high_resolution_clock::time_point TimeVar;
#define duration(a) std::chrono::duration_cast<std::chrono::microseconds>(a).count()
#define timeNow() std::chrono::high_resolution_clock::now()

//  We check that the values of the new regression match those of one built
//  from existing primitives.
TEST(ProbDistributionsBernoulliLogitGLM, glm_matches_bernoulli_logit_doubles) {
  Matrix<int, Dynamic, 1> n(3, 1);
  n << 1, 0, 1;
  Matrix<double, Dynamic, Dynamic> x(3, 2);
  x << -12, 46, -42,
      24, 25, 27;
  Matrix<double, Dynamic, 1> beta(2, 1);
  beta << 0.3, 2;
  double alpha = 0.3;
  Matrix<double, Dynamic, 1> alphavec = alpha * Matrix<double, 3, 1>::Ones();
  Matrix<double, Dynamic, 1> theta(3, 1);
  theta = x * beta + alphavec;

  EXPECT_FLOAT_EQ(
    (stan::math::bernoulli_logit_lpmf(n, theta)),
    (stan::math::bernoulli_logit_glm_lpmf(n, x, beta, alpha)));
  EXPECT_FLOAT_EQ(
    (stan::math::bernoulli_logit_lpmf<true>(n, theta)),
    (stan::math::bernoulli_logit_glm_lpmf<true>(n, x, beta, alpha)));
  EXPECT_FLOAT_EQ(
    (stan::math::bernoulli_logit_lpmf<false>(n, theta)),
    (stan::math::bernoulli_logit_glm_lpmf<false>(n, x, beta, alpha)));
  EXPECT_FLOAT_EQ(
    (stan::math::bernoulli_logit_lpmf<true, Matrix<int, Dynamic, 1>>
                                                        (n, theta)),
    (stan::math::bernoulli_logit_glm_lpmf<true, Matrix<int, Dynamic, 1>>
                                                        (n, x, beta, alpha)));
  EXPECT_FLOAT_EQ(
    (stan::math::bernoulli_logit_lpmf<false, Matrix<int, Dynamic, 1>>
                                                        (n, theta)),
    (stan::math::bernoulli_logit_glm_lpmf<false, Matrix<int, Dynamic, 1>>
                                                        (n, x, beta, alpha)));
  EXPECT_FLOAT_EQ(
    (stan::math::bernoulli_logit_lpmf<Matrix<int, Dynamic, 1>>
                                                        (n, theta)),
    (stan::math::bernoulli_logit_glm_lpmf<Matrix<int, Dynamic, 1>>
                                                        (n, x, beta, alpha)));
}

//  We check that the gradients of the new regression match those of one built
//  from existing primitives.
TEST(ProbDistributionsBernoulliLogitGLM, glm_matches_bernoulli_logit_vars) {
  Matrix<int, Dynamic, 1> n(3, 1);
  n << 1, 0, 1;
  Matrix<var, Dynamic, Dynamic> x(3, 2);
  x << -12, 46, -42,
      24, 25, 27;
  Matrix<var, Dynamic, 1> beta(2, 1);
  beta << 0.3, 2;
  var alpha = 0.3;
  Matrix<var, Dynamic, 1> alphavec = alpha * Matrix<double, 3, 1>::Ones();
  Matrix<var, Dynamic, 1> theta(3, 1);
  theta = x * beta + alphavec;

  var lp = stan::math::bernoulli_logit_lpmf(n, theta);
  lp.grad();

  stan::math::recover_memory();

  Matrix<int, Dynamic, 1> n2(3, 1);
  n2 << 1, 0, 1;
  Matrix<var, Dynamic, Dynamic> x2(3, 2);
  x2 << -12, 46, -42,
      24, 25, 27;
  Matrix<var, Dynamic, 1> beta2(2, 1);
  beta2 << 0.3, 2;
  var alpha2 = 0.3;

  var lp2 = stan::math::bernoulli_logit_glm_lpmf(n2, x2, beta2, alpha2);
  lp2.grad();

  EXPECT_FLOAT_EQ(lp.val(),
                  lp2.val());
  for (size_t i = 0; i < 2; i++) {
    EXPECT_FLOAT_EQ(beta[i].adj(), beta2[i].adj());
  }
  EXPECT_FLOAT_EQ(alpha.adj(), alpha2.adj());
  for (size_t j = 0; j < 3; j++) {
    for (size_t i = 0; i < 2; i++) {
      EXPECT_FLOAT_EQ(x(j, i).adj(), x2(j, i).adj());
    }
  }
}

//  Here, we compare the speed of the new regression to that of one built from
//  existing primitives.

/*
TEST(ProbDistributionsBernoulliLogitGLM, glm_matches_bernoulli_logit_speed) {
  const int R = 30000;
  const int C = 1000;

  Matrix<int, Dynamic, 1> n(R, 1);
  for (size_t i = 0; i < R; i++) {
    n[i] = rand()%2;
  }

  int T1 = 0;
  int T2 = 0;

  for (size_t testnumber = 0; testnumber < 30; testnumber++){
    Matrix<double, Dynamic, Dynamic> xreal = Matrix<double, Dynamic, Dynamic>::Random(R, C);
    Matrix<double, Dynamic, 1> betareal = Matrix<double, Dynamic, Dynamic>::Random(C, 1);
    Matrix<double, 1, 1> alphareal = Matrix<double, 1, 1>::Random(1, 1);
    Matrix<double, Dynamic, 1> alpharealvec = Matrix<double, R, 1>::Ones() * alphareal;

    Matrix<var, Dynamic, 1> beta = betareal;
    Matrix<var, Dynamic, 1> theta(R, 1);


    TimeVar t1 = timeNow();
    theta = (xreal * beta) + alpharealvec;
    var lp = stan::math::bernoulli_logit_lpmf(n, theta);

    lp.grad();
    TimeVar t2 = timeNow();

    stan::math::recover_memory();

    Matrix<var, Dynamic, 1> beta2 = betareal;

    TimeVar t3 = timeNow();
    var lp2 = stan::math::bernoulli_logit_glm_lpmf(n, xreal, beta2, alphareal[0]);
    lp2.grad();
    TimeVar t4 = timeNow();
    stan::math::recover_memory();
    T1 += duration(t2 - t1);
    T2 += duration(t4 - t3);

  }

  std::cout << "Existing Primitives:" << std::endl << T1 << std::endl  << "New Primitives:" << std::endl << T2 << std::endl;
}
*/
