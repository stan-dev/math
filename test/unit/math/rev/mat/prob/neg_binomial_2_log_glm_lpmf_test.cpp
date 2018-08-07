#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <stan/math/prim/scal/fun/value_of.hpp>

using Eigen::Dynamic;
using Eigen::Matrix;
using stan::math::var;

//  We check that the values of the new regression match those of one built
//  from existing primitives.
TEST(ProbDistributionsNegBinomial2LogGLM,
     glm_matches_neg_binomial_2_log_doubles) {
  Matrix<int, Dynamic, 1> n(3, 1);
  n << 1, 0, 1;
  Matrix<double, Dynamic, Dynamic> x(3, 2);
  x << -12, 46, -42, 24, 25, 27;
  Matrix<double, Dynamic, 1> beta(2, 1);
  beta << 0.3, 2;
  double alpha = 0.3;
  Matrix<double, Dynamic, 1> alphavec = alpha * Matrix<double, 3, 1>::Ones();
  Matrix<double, Dynamic, 1> theta(3, 1);
  theta = x * beta + alphavec;
  double phi = 2;
  EXPECT_FLOAT_EQ(
      (stan::math::neg_binomial_2_log_lpmf(n, theta, phi)),
      (stan::math::neg_binomial_2_log_glm_lpmf(n, x, alpha, beta, phi)));
  EXPECT_FLOAT_EQ(
      (stan::math::neg_binomial_2_log_lpmf<true>(n, theta, phi)),
      (stan::math::neg_binomial_2_log_glm_lpmf<true>(n, x, alpha, beta, phi)));
}
//  We check that the values of the new regression match those of one built
//  from existing primitives.
TEST(ProbDistributionsNegBinomial2LogGLM,
     glm_matches_neg_binomial_2_log_doubles_rand) {
  for (size_t ii = 0; ii < 20000; ii++) {
    Matrix<int, Dynamic, 1> n(3, 1);
    for (size_t i = 0; i < 3; i++) {
      n[i] = Matrix<uint, 1, 1>::Random(1, 1)[0] % 200;
    }
    Matrix<double, Dynamic, Dynamic> x
        = Matrix<double, Dynamic, Dynamic>::Random(3, 2);
    Matrix<double, Dynamic, 1> beta
        = Matrix<double, Dynamic, Dynamic>::Random(2, 1);
    Matrix<double, 1, 1> alphamat = Matrix<double, 1, 1>::Random(1, 1);
    double alpha = alphamat[0];
    Matrix<double, Dynamic, 1> alphavec = alpha * Matrix<double, 3, 1>::Ones();
    Matrix<double, Dynamic, 1> theta(3, 1);
    theta = x * beta + alphavec;
    double phi = Matrix<double, Dynamic, 1>::Random(1, 1)[0] + 1;
    EXPECT_FLOAT_EQ(
        (stan::math::neg_binomial_2_log_lpmf(n, theta, phi)),
        (stan::math::neg_binomial_2_log_glm_lpmf(n, x, alpha, beta, phi)));
    EXPECT_FLOAT_EQ((stan::math::neg_binomial_2_log_lpmf<true>(n, theta, phi)),
                    (stan::math::neg_binomial_2_log_glm_lpmf<true>(n, x, alpha,
                                                                   beta, phi)));
  }
}
//  We check that the gradients of the new regression match those of one built
//  from existing primitives.
TEST(ProbDistributionsNegBinomial2LogGLM, glm_matches_neg_binomial_2_log_vars) {
  Matrix<int, Dynamic, 1> n(3, 1);
  n << 1, 0, 1;
  Matrix<var, Dynamic, Dynamic> x(3, 2);
  x << -12, 46, -42, 24, 25, 27;
  Matrix<var, Dynamic, 1> beta(2, 1);
  beta << 0.3, 2;
  var alpha = 0.3;
  Matrix<var, Dynamic, 1> alphavec = alpha * Matrix<double, 3, 1>::Ones();
  Matrix<var, Dynamic, 1> theta(3, 1);
  theta = x * beta + alphavec;
  var phi = 2;
  var lp = stan::math::neg_binomial_2_log_lpmf(n, theta, phi);
  lp.grad();

  double lp_val = lp.val();
  double alpha_adj = alpha.adj();
  Matrix<double, Dynamic, Dynamic> x_adj(3, 2);
  Matrix<double, Dynamic, 1> beta_adj(2, 1);
  double phi_adj = phi.adj();
  for (size_t i = 0; i < 2; i++) {
    beta_adj[i] = beta[i].adj();
    for (size_t j = 0; j < 3; j++) {
      x_adj(j, i) = x(j, i).adj();
    }
  }

  stan::math::recover_memory();

  Matrix<int, Dynamic, 1> n2(3, 1);
  n2 << 1, 0, 1;
  Matrix<var, Dynamic, Dynamic> x2(3, 2);
  x2 << -12, 46, -42, 24, 25, 27;
  Matrix<var, Dynamic, 1> beta2(2, 1);
  beta2 << 0.3, 2;
  var alpha2 = 0.3;
  var phi2 = 2;
  var lp2
      = stan::math::neg_binomial_2_log_glm_lpmf(n2, x2, alpha2, beta2, phi2);
  lp2.grad();
  EXPECT_FLOAT_EQ(lp_val, lp2.val());
  for (size_t i = 0; i < 2; i++) {
    EXPECT_FLOAT_EQ(beta_adj[i], beta2[i].adj());
  }
  EXPECT_FLOAT_EQ(alpha_adj, alpha2.adj());
  EXPECT_FLOAT_EQ(phi_adj, phi2.adj());
  for (size_t j = 0; j < 3; j++) {
    for (size_t i = 0; i < 2; i++) {
      EXPECT_FLOAT_EQ(x_adj(j, i), x2(j, i).adj());
    }
  }
}

//  We check that the gradients of the new regression match those of one built
//  from existing primitives.
TEST(ProbDistributionsNegBinomial2LogGLM,
     glm_matches_neg_binomial_2_log_vars_rand) {
  for (size_t ii = 0; ii < 200; ii++) {
    Matrix<int, Dynamic, 1> n(3, 1);
    for (size_t i = 0; i < 3; i++) {
      n[i] = Matrix<uint, 1, 1>::Random(1, 1)[0] % 200;
    }
    Matrix<double, Dynamic, Dynamic> xreal
        = Matrix<double, Dynamic, Dynamic>::Random(3, 2);
    Matrix<double, Dynamic, 1> betareal
        = Matrix<double, Dynamic, Dynamic>::Random(2, 1);
    Matrix<double, 1, 1> alphareal = Matrix<double, 1, 1>::Random(1, 1);
    double phireal = Matrix<double, Dynamic, 1>::Random(1, 1)[0] + 1;
    Matrix<var, Dynamic, 1> beta = betareal;
    Matrix<var, Dynamic, 1> theta(3, 1);
    Matrix<var, Dynamic, Dynamic> x = xreal;
    var alpha = alphareal[0];
    Matrix<var, Dynamic, 1> alphavec = Matrix<double, 3, 1>::Ones() * alpha;
    var phi = phireal;
    theta = (x * beta) + alphavec;
    var lp = stan::math::neg_binomial_2_log_lpmf(n, theta, phi);
    lp.grad();

    double lp_val = lp.val();
    double alpha_adj = alpha.adj();
    Matrix<double, Dynamic, Dynamic> x_adj(3, 2);
    Matrix<double, Dynamic, 1> beta_adj(2, 1);
    for (size_t i = 0; i < 2; i++) {
      beta_adj[i] = beta[i].adj();
      for (size_t j = 0; j < 3; j++) {
        x_adj(j, i) = x(j, i).adj();
      }
    }
    double phi_adj = phi.adj();

    stan::math::recover_memory();
    Matrix<var, Dynamic, 1> beta2 = betareal;
    Matrix<var, Dynamic, Dynamic> x2 = xreal;
    var alpha2 = alphareal[0];
    var phi2 = phireal;
    var lp2
        = stan::math::neg_binomial_2_log_glm_lpmf(n, x2, alpha2, beta2, phi2);
    lp2.grad();
    EXPECT_FLOAT_EQ(lp_val, lp2.val());
    for (size_t i = 0; i < 2; i++) {
      EXPECT_FLOAT_EQ(beta_adj[i], beta2[i].adj());
    }
    EXPECT_FLOAT_EQ(alpha_adj, alpha2.adj());
    EXPECT_FLOAT_EQ(phi_adj, phi2.adj());
    for (size_t j = 0; j < 3; j++) {
      for (size_t i = 0; i < 2; i++) {
        EXPECT_FLOAT_EQ(x_adj(j, i), x2(j, i).adj());
      }
    }
  }
}

//  We check that the gradients of the new regression match those of one built
//  from existing primitives.
TEST(ProbDistributionsNegBinomial2LogGLM,
     glm_matches_neg_binomial_2_log_varying_intercept) {
  for (size_t ii = 0; ii < 200; ii++) {
    Matrix<int, Dynamic, 1> n(3, 1);
    for (size_t i = 0; i < 3; i++) {
      n[i] = Matrix<uint, 1, 1>::Random(1, 1)[0] % 200;
    }
    Matrix<double, Dynamic, Dynamic> xreal
        = Matrix<double, Dynamic, Dynamic>::Random(3, 2);
    Matrix<double, Dynamic, 1> betareal
        = Matrix<double, Dynamic, Dynamic>::Random(2, 1);
    Matrix<double, Dynamic, 1> alphareal
        = Matrix<double, Dynamic, Dynamic>::Random(3, 1);
    Matrix<double, Dynamic, 1> phireal
        = Matrix<double, Dynamic, Dynamic>::Random(3, 1)
          + Matrix<double, Dynamic, 1>::Ones(3, 1);
    Matrix<var, Dynamic, 1> beta = betareal;
    Matrix<var, Dynamic, 1> theta(3, 1);
    Matrix<var, Dynamic, Dynamic> x = xreal;
    Matrix<var, Dynamic, 1> alpha = alphareal;
    Matrix<var, Dynamic, 1> phi = phireal;
    theta = (x * beta) + alpha;
    var lp = stan::math::neg_binomial_2_log_lpmf(n, theta, phi);
    lp.grad();

    double lp_val = lp.val();
    Matrix<double, Dynamic, 1> alpha_adj(3, 1);
    Matrix<double, Dynamic, Dynamic> x_adj(3, 2);
    Matrix<double, Dynamic, 1> beta_adj(2, 1);
    Matrix<double, Dynamic, 1> phi_adj(3, 1);
    for (size_t i = 0; i < 2; i++) {
      beta_adj[i] = beta[i].adj();
      for (size_t j = 0; j < 3; j++) {
        x_adj(j, i) = x(j, i).adj();
      }
    }
    for (size_t j = 0; j < 3; j++) {
      phi_adj[j] = phi[j].adj();
      alpha_adj[j] = alpha[j].adj();
    }

    stan::math::recover_memory();
    Matrix<var, Dynamic, 1> beta2 = betareal;
    Matrix<var, Dynamic, Dynamic> x2 = xreal;
    Matrix<var, Dynamic, 1> alpha2 = alphareal;
    Matrix<var, Dynamic, 1> phi2 = phireal;
    var lp2
        = stan::math::neg_binomial_2_log_glm_lpmf(n, x2, alpha2, beta2, phi2);
    lp2.grad();
    EXPECT_FLOAT_EQ(lp_val, lp2.val());
    for (size_t i = 0; i < 2; i++) {
      EXPECT_FLOAT_EQ(beta_adj[i], beta2[i].adj());
    }
    for (size_t j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(phi_adj[j], phi2[j].adj());
      EXPECT_FLOAT_EQ(alpha_adj[j], alpha2[j].adj());
      for (size_t i = 0; i < 2; i++) {
        EXPECT_FLOAT_EQ(x_adj(j, i), x2(j, i).adj());
      }
    }
  }
}
