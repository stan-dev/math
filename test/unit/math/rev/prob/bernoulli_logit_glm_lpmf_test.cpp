#include <stan/math/rev.hpp>
#include <stan/math/prim.hpp>
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>
#include <vector>
#include <cmath>

//  We check that the values of the new regression match those of one built
//  from existing primitives.
TEST(ProbDistributionsBernoulliLogitGLM, glm_matches_bernoulli_logit_doubles) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::var;
  using std::vector;
  vector<int> y{1, 0, 1};
  Matrix<double, Dynamic, Dynamic> x(3, 2);
  x << -12, 46, -42, 24, 25, 27;
  Matrix<double, Dynamic, 1> beta(2, 1);
  beta << 0.3, 2;
  double alpha = 0.3;
  Matrix<double, Dynamic, 1> alphavec = alpha * Matrix<double, 3, 1>::Ones();
  Matrix<double, Dynamic, 1> theta(3, 1);
  theta = x * beta + alphavec;

  EXPECT_FLOAT_EQ((stan::math::bernoulli_logit_lpmf(y, theta)),
                  (stan::math::bernoulli_logit_glm_lpmf(y, x, alpha, beta)));
  EXPECT_FLOAT_EQ(
      (stan::math::bernoulli_logit_lpmf<true>(y, theta)),
      (stan::math::bernoulli_logit_glm_lpmf<true>(y, x, alpha, beta)));
}

//  We check that the values of the new regression match those of one built
//  from existing primitives.
TEST(ProbDistributionsBernoulliLogitGLM,
     glm_matches_bernoulli_logit_doubles_rand) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::var;
  using std::vector;
  for (size_t ii = 0; ii < 200; ii++) {
    vector<int> y(3);
    for (size_t i = 0; i < 3; i++) {
      y[i] = Matrix<unsigned int, 1, 1>::Random(1, 1)[0] % 2;
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

    EXPECT_FLOAT_EQ((stan::math::bernoulli_logit_lpmf(y, theta)),
                    (stan::math::bernoulli_logit_glm_lpmf(y, x, alpha, beta)));
    EXPECT_FLOAT_EQ(
        (stan::math::bernoulli_logit_lpmf<true>(y, theta)),
        (stan::math::bernoulli_logit_glm_lpmf<true>(y, x, alpha, beta)));
  }
}

namespace stan {
namespace test {
namespace bernoulli_glm {
//  We check that the gradients of the new regression match those of one built
//  from existing primitives.
template <bool UseVarMat>
void matches_bernoulli_logit_vars() {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::var;
  using std::vector;
  vector<int> y{1, 0, 1};
  Matrix<var, Dynamic, Dynamic> x(3, 2);
  x << -1234, 46, -42, 24, 25, 27;
  Matrix<var, Dynamic, 1> beta(2, 1);
  beta << 0.3, 2000;
  var alpha = 0.3;
  Matrix<var, Dynamic, 1> alphavec = alpha * Matrix<double, 3, 1>::Ones();
  Matrix<var, Dynamic, 1> theta(3, 1);
  theta = x * beta + alphavec;

  var lp = stan::math::bernoulli_logit_lpmf(y, theta);
  lp.grad();

  double lp_val = lp.val();
  double alpha_adj = alpha.adj();
  Matrix<double, Dynamic, Dynamic> x_adj(3, 2);
  Matrix<double, Dynamic, 1> beta_adj(2, 1);
  for (size_t i = 0; i < 2; i++) {
    beta_adj[i] = beta.adj()[i];
    for (size_t j = 0; j < 3; j++) {
      x_adj(j, i) = x.adj()(j, i);
    }
  }

  stan::math::recover_memory();
  using stan::math::test::internal::cond_matrix;
  vector<int> y2{1, 0, 1};
  Matrix<double, Dynamic, Dynamic> x2_val(3, 2);
  x2_val << -1234, 46, -42, 24, 25, 27;
  cond_matrix<UseVarMat, Matrix<double, Dynamic, Dynamic>> x2 = x2_val;
  Matrix<double, Dynamic, 1> beta2_val(2, 1);
  beta2_val << 0.3, 2000;
  cond_matrix<UseVarMat, Matrix<double, Dynamic, 1>> beta2 = beta2_val;
  var alpha2 = 0.3;

  var lp2 = stan::math::bernoulli_logit_glm_lpmf(y2, x2, alpha2, beta2);
  lp2.grad();

  EXPECT_FLOAT_EQ(lp_val, lp2.val());
  EXPECT_FLOAT_EQ(alpha_adj, alpha2.adj());
  for (size_t i = 0; i < 2; i++) {
    EXPECT_FLOAT_EQ(beta_adj[i], beta2.adj().coeff(i));
    for (size_t j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(x_adj(j, i), x2.adj().coeff(j, i));
    }
  }
}


template <bool UseVarMat>
void broadcast_x() {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::var;
  using std::vector;
  using stan::math::test::internal::cond_matrix;
  vector<int> y{1, 0, 1};
  Matrix<double, 1, Dynamic> x(1, 2);
  x << -12, 46;
  cond_matrix<UseVarMat, Matrix<double, 1, Dynamic>> x1 = x;
  cond_matrix<UseVarMat, Matrix<double, Dynamic, Dynamic>> x_mat = x.replicate(3, 1);
  Matrix<double, Dynamic, 1> beta(2, 1);
  beta << 0.3, 2;
  cond_matrix<UseVarMat, Matrix<double, Dynamic, 1>> beta1 = beta;
  cond_matrix<UseVarMat, Matrix<double, Dynamic, 1>> beta2 = beta;
  var alpha1 = 0.3;
  var alpha2 = 0.3;

  var lp1 = stan::math::bernoulli_logit_glm_lpmf(y, x1, alpha1, beta1);
  var lp2 = stan::math::bernoulli_logit_glm_lpmf(y, x_mat, alpha2, beta2);

  EXPECT_DOUBLE_EQ(lp1.val(), lp2.val());

  (lp1 + lp2).grad();

  for (int i = 0; i < 2; i++) {
    EXPECT_DOUBLE_EQ(x1.adj()[i], x_mat.col(i).adj().sum());
    EXPECT_DOUBLE_EQ(beta1.adj()[i], beta2.adj()[i]);
  }
  EXPECT_DOUBLE_EQ(alpha1.adj(), alpha2.adj());
}


template <bool UseVarMat>
void broadcast_y() {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::var;
  using std::vector;
  using stan::math::test::internal::cond_matrix;
  int y = 1;
  Matrix<int, Dynamic, 1> y_vec = Matrix<int, Dynamic, 1>::Constant(3, 1, y);
  Matrix<double, Dynamic, Dynamic> x(3, 2);
  x << -12, 46, -42, 24, 25, 27;
  cond_matrix<UseVarMat, Matrix<double, Dynamic, Dynamic>> x1 = x;
  cond_matrix<UseVarMat, Matrix<double, Dynamic, Dynamic>> x2 = x;
  Matrix<double, Dynamic, 1> beta(2, 1);
  beta << 0.3, 2;
  cond_matrix<UseVarMat, Matrix<double, Dynamic, 1>> beta1 = beta;
  cond_matrix<UseVarMat, Matrix<double, Dynamic, 1>> beta2 = beta;
  var alpha1 = 0.3;
  var alpha2 = 0.3;

  var lp1 = stan::math::bernoulli_logit_glm_lpmf(y, x1, alpha1, beta1);
  var lp2 = stan::math::bernoulli_logit_glm_lpmf(y_vec, x2, alpha2, beta2);

  EXPECT_DOUBLE_EQ(lp1.val(), lp2.val());

  (lp1 + lp2).grad();

  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      EXPECT_DOUBLE_EQ(x1.adj()(j, i), x2.adj()(j, i));
    }
    EXPECT_DOUBLE_EQ(beta1.adj()[i], beta2.adj()[i]);
  }
  EXPECT_DOUBLE_EQ(alpha1.adj(), alpha2.adj());
}


template <bool UseVarMat>
void matches_bernoulli_logit_vars_zero_instances() {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::var;
  using std::vector;
  using stan::math::test::internal::cond_matrix;
  vector<int> y{};
  Matrix<var, Dynamic, Dynamic> x(0, 2);
  Matrix<var, Dynamic, 1> beta(2, 1);
  beta << 0.3, 2000;
  var alpha = 0.3;
  Matrix<var, Dynamic, 1> theta(0, 1);

  var lp = stan::math::bernoulli_logit_lpmf(y, theta);
  lp.grad();

  double lp_val = lp.val();
  double alpha_adj = alpha.adj();
  Matrix<double, Dynamic, 1> beta_adj(2, 1);
  for (size_t i = 0; i < 2; i++) {
    beta_adj[i] = beta.adj()[i];
  }

  stan::math::recover_memory();

  vector<int> y2{};
  Matrix<double, Dynamic, Dynamic> x2_val(0, 2);
  Matrix<double, Dynamic, 1> beta2_val(2, 1);
  beta2_val << 0.3, 2000;
  var alpha2 = 0.3;
  cond_matrix<UseVarMat, Matrix<double, Dynamic, Dynamic>> x2 = x2_val;
  cond_matrix<UseVarMat, Matrix<double, Dynamic, 1>> beta2 = beta2_val;
  var lp2 = stan::math::bernoulli_logit_glm_lpmf(y2, x2, alpha2, beta2);
  lp2.grad();

  EXPECT_FLOAT_EQ(lp_val, lp2.val());
  EXPECT_FLOAT_EQ(alpha_adj, alpha2.adj());
  for (size_t i = 0; i < 2; i++) {
    EXPECT_FLOAT_EQ(beta_adj[i], beta2.adj()[i]);
  }
}


template <bool UseVarMat>
void matches_bernoulli_logit_vars_zero_attributes() {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::var;
  using std::vector;
  using stan::math::test::internal::cond_matrix;
  vector<int> y{1, 0, 1};
  Matrix<var, Dynamic, Dynamic> x(3, 0);
  Matrix<var, Dynamic, 1> beta(0, 1);
  var alpha = 0.3;
  Matrix<var, Dynamic, 1> alphavec = alpha * Matrix<double, 3, 1>::Ones();
  Matrix<var, Dynamic, 1> theta(3, 1);
  theta = x * beta + alphavec;

  var lp = stan::math::bernoulli_logit_lpmf(y, theta);
  lp.grad();

  double lp_val = lp.val();
  double alpha_adj = alpha.adj();

  stan::math::recover_memory();

  vector<int> y2{1, 0, 1};
  Matrix<var, Dynamic, Dynamic> x2(3, 0);
  Matrix<var, Dynamic, 1> beta2(0, 1);
  var alpha2 = 0.3;

  var lp2 = stan::math::bernoulli_logit_glm_lpmf(y2, x2, alpha2, beta2);
  lp2.grad();

  EXPECT_FLOAT_EQ(lp_val, lp2.val());
  EXPECT_FLOAT_EQ(alpha_adj, alpha2.adj());
}


//  We check that the gradients of the new regression match those of one built
//  from existing primitives.
template <bool UseVarMat>
void matches_bernoulli_logit_vars_rand() {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::var;
  using std::vector;
  using stan::math::test::internal::cond_matrix;
  for (size_t ii = 0; ii < 42; ii++) {
    vector<int> y(3);
    for (size_t i = 0; i < 3; i++) {
      y[i] = Matrix<unsigned int, 1, 1>::Random(1, 1)[0] % 2;
    }
    Matrix<double, Dynamic, Dynamic> xreal
        = Matrix<double, Dynamic, Dynamic>::Random(3, 2);
    Matrix<double, Dynamic, 1> betareal
        = Matrix<double, Dynamic, Dynamic>::Random(2, 1);
    Matrix<double, 1, 1> alphareal = Matrix<double, 1, 1>::Random(1, 1);

    Matrix<var, Dynamic, 1> beta = betareal;
    Matrix<var, Dynamic, 1> theta(3, 1);
    Matrix<var, Dynamic, Dynamic> x = xreal;
    var alpha = alphareal[0];
    Matrix<var, Dynamic, 1> alphavec = Matrix<double, 3, 1>::Ones() * alpha;

    theta = (x * beta) + alphavec;
    var lp = stan::math::bernoulli_logit_lpmf(y, theta);

    lp.grad();

    double lp_val = lp.val();
    double alpha_adj = alpha.adj();
    Matrix<double, Dynamic, Dynamic> x_adj(3, 2);
    Matrix<double, Dynamic, 1> beta_adj(2, 1);
    for (size_t i = 0; i < 2; i++) {
      beta_adj[i] = beta.adj()[i];
      for (size_t j = 0; j < 3; j++) {
        x_adj(j, i) = x.adj()(j, i);
      }
    }

    stan::math::recover_memory();

    Matrix<var, Dynamic, 1> beta2 = betareal;
    Matrix<var, Dynamic, Dynamic> x2 = xreal;
    var alpha2 = alphareal[0];

    var lp2 = stan::math::bernoulli_logit_glm_lpmf(y, x2, alpha2, beta2);
    lp2.grad();

    EXPECT_FLOAT_EQ(lp_val, lp2.val());
    EXPECT_FLOAT_EQ(alpha_adj, alpha2.adj());
    for (size_t i = 0; i < 2; i++) {
      EXPECT_FLOAT_EQ(beta_adj[i], beta2.adj()[i]);
      for (size_t j = 0; j < 3; j++) {
        EXPECT_FLOAT_EQ(x_adj(j, i), x2.adj()(j, i));
      }
    }
  }
}


//  We check that the gradients of the new regression match those of one built
//  from existing primitives, in case beta is a scalar.
template <bool UseVarMat>
void matches_bernoulli_logit_vars_rand_scal_beta() {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::var;
  using std::vector;
  using stan::math::test::internal::cond_matrix;
  for (size_t ii = 0; ii < 42; ii++) {
    vector<int> y(3);
    for (size_t i = 0; i < 3; i++) {
      y[i] = Matrix<unsigned int, 1, 1>::Random(1, 1)[0] % 2;
    }
    Matrix<double, Dynamic, Dynamic> xreal
        = Matrix<double, Dynamic, Dynamic>::Random(3, 1);
    double betareal = Matrix<double, Dynamic, Dynamic>::Random(1, 1)(0, 0);
    Matrix<double, 1, 1> alphareal = Matrix<double, 1, 1>::Random(1, 1);
    var beta = betareal;
    Matrix<var, Dynamic, 1> theta(3, 1);
    Matrix<var, Dynamic, Dynamic> x = xreal;
    var alpha = alphareal[0];
    Matrix<var, Dynamic, 1> alphavec = Matrix<double, 3, 1>::Ones() * alpha;
    theta = (x * beta) + alphavec;
    var lp = stan::math::bernoulli_logit_lpmf(y, theta);
    lp.grad();

    double lp_val = lp.val();
    double alpha_adj = alpha.adj();
    Matrix<double, Dynamic, Dynamic> x_adj(3, 1);
    double beta_adj = beta.adj();
    for (size_t j = 0; j < 3; j++) {
      x_adj(j, 0) = x(j, 0).adj();
    }

    stan::math::recover_memory();

    var beta2 = betareal;
    Matrix<var, Dynamic, Dynamic> x2 = xreal;
    var alpha2 = alphareal[0];
    var lp2 = stan::math::bernoulli_logit_glm_lpmf(y, x2, alpha2, beta2);
    lp2.grad();
    EXPECT_FLOAT_EQ(lp_val, lp2.val());
    EXPECT_FLOAT_EQ(beta_adj, beta2.adj());
    EXPECT_FLOAT_EQ(alpha_adj, alpha2.adj());
    for (size_t j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(x_adj(j, 0), x2(j, 0).adj());
    }
  }
}


//  We check that the gradients of the new regression match those of one built
//  from existing primitives, for the GLM with varying intercept.
template <bool UseVarMat>
void matches_bernoulli_varying_intercept() {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::var;
  using std::vector;
  using stan::math::test::internal::cond_matrix;
  for (size_t ii = 0; ii < 42; ii++) {
    vector<int> y(3);
    for (size_t i = 0; i < 3; i++) {
      y[i] = Matrix<unsigned int, 1, 1>::Random(1, 1)[0] % 2;
    }
    Matrix<double, Dynamic, Dynamic> xreal
        = Matrix<double, Dynamic, Dynamic>::Random(3, 2);
    Matrix<double, Dynamic, 1> betareal
        = Matrix<double, Dynamic, Dynamic>::Random(2, 1);
    Matrix<double, Dynamic, 1> alphareal
        = Matrix<double, Dynamic, Dynamic>::Random(3, 1);

    Matrix<var, Dynamic, 1> beta = betareal;
    Matrix<var, Dynamic, 1> theta(3, 1);
    Matrix<var, Dynamic, Dynamic> x = xreal;
    Matrix<var, Dynamic, 1> alpha = alphareal;
    theta = (x * beta) + alpha;
    var lp = stan::math::bernoulli_logit_lpmf(y, theta);

    lp.grad();

    double lp_val = lp.val();
    Matrix<double, Dynamic, 1> alpha_adj(3, 1);
    Matrix<double, Dynamic, Dynamic> x_adj(3, 2);
    Matrix<double, Dynamic, 1> beta_adj(2, 1);
    for (size_t i = 0; i < 2; i++) {
      beta_adj[i] = beta.adj()[i];
      for (size_t j = 0; j < 3; j++) {
        x_adj(j, i) = x.adj()(j, i);
      }
    }
    for (size_t j = 0; j < 3; j++) {
      alpha_adj[j] = alpha[j].adj();
    }

    stan::math::recover_memory();

    Matrix<var, Dynamic, 1> beta2 = betareal;
    Matrix<var, Dynamic, Dynamic> x2 = xreal;
    Matrix<var, Dynamic, 1> alpha2 = alphareal;

    var lp2 = stan::math::bernoulli_logit_glm_lpmf(y, x2, alpha2, beta2);
    lp2.grad();

    EXPECT_FLOAT_EQ(lp_val, lp2.val());
    for (size_t i = 0; i < 2; i++) {
      EXPECT_FLOAT_EQ(beta_adj[i], beta2.adj()[i]);
      for (size_t j = 0; j < 3; j++) {
        EXPECT_FLOAT_EQ(x_adj(j, i), x2.adj()(j, i));
      }
    }
    for (size_t j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(alpha_adj[j], alpha2[j].adj());
    }
  }
}

//  We check that we can instantiate all different interface types.
template <bool UseVarMat>
void matches_bernoulli_logit_interface_types() {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::var;
  using std::vector;
  using stan::math::test::internal::cond_matrix;
  double value = 0;
  double value2 = 0;

  int i = 1;
  std::vector<int> vi = {{1, 0}};
  double d = 1.0;
  std::vector<double> vd = {{1.0, 2.0}};
  Eigen::VectorXd ev(2);
  Eigen::RowVectorXd rv(2);
  Eigen::MatrixXd m1(1, 1);
  m1 << 1.0;
  Eigen::MatrixXd m(2, 2);
  ev << 1.0, 2.0;
  rv << 1.0, 2.0;
  m << 1.0, 2.0, 3.0, 4.0;

  value += stan::math::bernoulli_logit_glm_lpmf(i, m1, d, d);
  value += stan::math::bernoulli_logit_glm_lpmf(vi, m, vd, vd);
  value += stan::math::bernoulli_logit_glm_lpmf(vi, m, ev, ev);
  value += stan::math::bernoulli_logit_glm_lpmf(vi, m, rv, rv);

  var v = 1.0;
  std::vector<var> vv = {{1.0, 2.0}};
  Eigen::Matrix<double, -1, 1> evv_val(2);
  evv_val << 1.0, 2.0;
  cond_matrix<UseVarMat, Eigen::Matrix<double, -1, 1>> evv = evv_val;
  Eigen::Matrix<double, 1, -1> rvv_val(2);
  rvv_val << 1.0, 2.0;
  cond_matrix<UseVarMat, Eigen::Matrix<double, 1, -1>> rvv = rvv_val;
  Eigen::Matrix<double, -1, -1> m1v_val(1, 1);
  m1v_val << 1.0;
  cond_matrix<UseVarMat, Eigen::Matrix<double, -1, -1>> m1v = m1v_val;
  Eigen::Matrix<double, -1, -1> mv_val(2, 2);
  mv_val << 1.0, 2.0, 3.0, 4.0;
  cond_matrix<UseVarMat, Eigen::Matrix<double, -1, -1>> mv = mv_val;

  value2 += stan::math::bernoulli_logit_glm_lpmf(i, m1v, v, v).val();
  value2 += stan::math::bernoulli_logit_glm_lpmf(vi, mv, vv, vv).val();
  value2 += stan::math::bernoulli_logit_glm_lpmf(vi, mv, evv, evv).val();
  value2 += stan::math::bernoulli_logit_glm_lpmf(vi, mv, rvv, rvv).val();

  EXPECT_FLOAT_EQ(value, value2);
}

}
}
}

TEST(ProbDistributionsBernoulliLogitGLM, glm_matches_bernoulli_logit_vars) {
 stan::test::bernoulli_glm::matches_bernoulli_logit_vars<true>();
 stan::test::bernoulli_glm::matches_bernoulli_logit_vars<false>();
}

TEST(ProbDistributionsBernoulliLogitGLM, broadcast_x) {
  stan::test::bernoulli_glm::broadcast_x<true>();
  stan::test::bernoulli_glm::broadcast_x<false>();
}

TEST(ProbDistributionsBernoulliLogitGLM, broadcast_y) {
  stan::test::bernoulli_glm::broadcast_y<true>();
  stan::test::bernoulli_glm::broadcast_y<false>();
}

TEST(ProbDistributionsBernoulliLogitGLM,
     glm_matches_bernoulli_logit_vars_zero_instances) {
stan::test::bernoulli_glm::matches_bernoulli_logit_vars_zero_instances<true>();
stan::test::bernoulli_glm::matches_bernoulli_logit_vars_zero_instances<false>();
}

TEST(ProbDistributionsBernoulliLogitGLM,
     glm_matches_bernoulli_logit_vars_zero_attributes) {
stan::test::bernoulli_glm::matches_bernoulli_logit_vars_zero_attributes<true>();
stan::test::bernoulli_glm::matches_bernoulli_logit_vars_zero_attributes<false>();
}

TEST(ProbDistributionsBernoulliLogitGLM,
     glm_matches_bernoulli_logit_vars_rand) {
stan::test::bernoulli_glm::matches_bernoulli_logit_vars_rand<true>();
stan::test::bernoulli_glm::matches_bernoulli_logit_vars_rand<false>();
}

TEST(ProbDistributionsBernoulliLogitGLM,
     glm_matches_bernoulli_logit_vars_rand_scal_beta) {
stan::test::bernoulli_glm::matches_bernoulli_logit_vars_rand_scal_beta<true>();
stan::test::bernoulli_glm::matches_bernoulli_logit_vars_rand_scal_beta<false>();
}


TEST(ProbDistributionsBernoulliLogitGLM,
     glm_matches_bernoulli_varying_intercept) {
stan::test::bernoulli_glm::matches_bernoulli_varying_intercept<true>();
stan::test::bernoulli_glm::matches_bernoulli_varying_intercept<false>();
}

TEST(ProbDistributionsBernoulliLogitGLM,
     glm_matches_bernoulli_logit_interface_types) {
stan::test::bernoulli_glm::matches_bernoulli_logit_interface_types<true>();
stan::test::bernoulli_glm::matches_bernoulli_logit_interface_types<false>();
}

//  We check that the right errors are thrown.
TEST(ProbDistributionsBernoulliLogitGLM,
     glm_matches_bernoulli_logit_error_checking) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::var;
  using std::vector;
  int N = 3;
  int M = 2;
  int W = 4;

  Eigen::Matrix<int, -1, 1> y(N, 1);
  for (int n = 0; n < N; n++) {
    y[n] = Eigen::Matrix<unsigned int, -1, 1>::Random(1, 1)[0] % 2;
  }
  Eigen::Matrix<int, -1, 1> yw1(W, 1);
  for (int n = 0; n < W; n++) {
    yw1[n] = Eigen::Matrix<unsigned int, -1, 1>::Random(1, 1)[0] % 2;
  }
  Eigen::Matrix<int, -1, 1> yw2(N, 1);
  for (int n = 0; n < N; n++) {
    yw2[n] = -(Eigen::Matrix<unsigned int, -1, 1>::Random(1, 1)[0] % 2);
  }
  Eigen::Matrix<int, -1, 1> yw3(N, 1);
  for (int n = 0; n < N; n++) {
    yw3[n] = 42 + (Eigen::Matrix<unsigned int, -1, 1>::Random(1, 1)[0] % 2);
  }
  Eigen::Matrix<double, -1, -1> x = Eigen::Matrix<double, -1, -1>::Random(N, M);
  Eigen::Matrix<double, -1, -1> xw1
      = Eigen::Matrix<double, -1, -1>::Random(W, M);
  Eigen::Matrix<double, -1, -1> xw2
      = Eigen::Matrix<double, -1, -1>::Random(N, W);
  Eigen::Matrix<double, -1, -1> xw3
      = Eigen::Matrix<double, -1, -1>::Random(N, M) * NAN;
  Eigen::Matrix<double, -1, 1> alpha
      = Eigen::Matrix<double, -1, 1>::Random(N, 1);
  Eigen::Matrix<double, -1, 1> alphaw1
      = Eigen::Matrix<double, -1, 1>::Random(W, 1);
  Eigen::Matrix<double, -1, 1> alphaw2
      = Eigen::Matrix<double, -1, 1>::Random(N, 1) * NAN;
  Eigen::Matrix<double, -1, 1> beta
      = Eigen::Matrix<double, -1, 1>::Random(M, 1);
  Eigen::Matrix<double, -1, 1> betaw1
      = Eigen::Matrix<double, -1, 1>::Random(W, 1);
  Eigen::Matrix<double, -1, 1> betaw2
      = Eigen::Matrix<double, -1, 1>::Random(M, 1) * NAN;

  EXPECT_THROW(stan::math::bernoulli_logit_glm_lpmf(yw1, x, alpha, beta),
               std::invalid_argument);
  EXPECT_THROW(stan::math::bernoulli_logit_glm_lpmf(yw2, x, alpha, beta),
               std::domain_error);
  EXPECT_THROW(stan::math::bernoulli_logit_glm_lpmf(yw3, x, alpha, beta),
               std::domain_error);
  EXPECT_THROW(stan::math::bernoulli_logit_glm_lpmf(y, xw1, alpha, beta),
               std::invalid_argument);
  EXPECT_THROW(stan::math::bernoulli_logit_glm_lpmf(y, xw2, alpha, beta),
               std::invalid_argument);
  EXPECT_THROW(stan::math::bernoulli_logit_glm_lpmf(y, xw3, alpha, beta),
               std::domain_error);
  EXPECT_THROW(stan::math::bernoulli_logit_glm_lpmf(y, x, alphaw1, beta),
               std::invalid_argument);
  EXPECT_THROW(stan::math::bernoulli_logit_glm_lpmf(y, x, alphaw2, beta),
               std::domain_error);
  EXPECT_THROW(stan::math::bernoulli_logit_glm_lpmf(y, x, alpha, betaw1),
               std::invalid_argument);
  EXPECT_THROW(stan::math::bernoulli_logit_glm_lpmf(y, x, alpha, betaw2),
               std::domain_error);
}
