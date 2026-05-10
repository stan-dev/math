#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <vector>

// -----------------------------------------------------------------------
// Equivalence with multinomial_logit_lpmf (single instance, N=1)
// -----------------------------------------------------------------------
TEST(ProbMultinomialLogitGLM, matchesMultinomialLogit) {
  // N=1, K=4 classes, M=3 features
  Eigen::MatrixXd x(1, 3);
  x << 0.5, -1.2, 0.8;
  Eigen::RowVectorXd alpha(4);
  alpha << 0.1, -0.3, 0.5, -0.2;
  Eigen::MatrixXd beta(3, 4);
  beta <<  0.2, -0.1,  0.4,  0.0,
          -0.3,  0.5, -0.2,  0.1,
           0.1,  0.0,  0.3, -0.4;

  // eta = alpha^T + beta^T x^T  (K-vector for multinomial_logit_lpmf)
  Eigen::VectorXd eta = alpha.transpose() + beta.transpose() * x.transpose();

  std::vector<int> ns{3, 1, 0, 2};
  std::vector<std::vector<int>> y{ns};

  EXPECT_FLOAT_EQ(stan::math::multinomial_logit_lpmf(ns, eta),
                  stan::math::multinomial_logit_glm_lpmf(y, x, alpha, beta));
}

// -----------------------------------------------------------------------
// Equivalence with binomial_logit_lpmf (K=2 special case)
//
// With K=2, alpha=[0,0], and beta_mult = [beta_bin | 0_col] (M x 2),
// the linear predictor per instance is lin_n = [x_n·beta_bin, 0].
// The softmax of [a, 0] is [sigmoid(a), sigmoid(-a)], so the
// multinomial log-PMF reduces to binomial_logit_lpmf(y0 | y0+y1, x·beta_bin).
// -----------------------------------------------------------------------
TEST(ProbMultinomialLogitGLM, matchesBinomialLogitLpmf) {
  // N=2 instances, M=2 features, K=2 classes
  Eigen::MatrixXd x(2, 2);
  x << 1.0, -0.5, 0.3, 0.7;

  Eigen::VectorXd beta_bin(2);
  beta_bin << 0.4, -0.6;

  // Multinomial parameters: alpha=[0,0], beta_mult = [[beta_bin], [0]]
  const Eigen::RowVectorXd alpha_mult = Eigen::RowVectorXd::Zero(2);
  Eigen::MatrixXd beta_mult(2, 2);
  beta_mult << beta_bin(0), 0.0,
               beta_bin(1), 0.0;

  const std::vector<int> n_success{3, 2};
  const std::vector<int> n_trials{5, 4};

  std::vector<std::vector<int>> y(2);
  for (int i = 0; i < 2; ++i) {
    y[i] = {n_success[i], n_trials[i] - n_success[i]};
  }

  // Binomial logit uses eta = x * beta_bin (no intercept, matches alpha=[0,0])
  const Eigen::VectorXd eta = x * beta_bin;
  EXPECT_FLOAT_EQ(stan::math::binomial_logit_lpmf(n_success, n_trials, eta),
                  stan::math::multinomial_logit_glm_lpmf(y, x, alpha_mult,
                                                         beta_mult));
}

// -----------------------------------------------------------------------
// K=3, N=3 instances, M=2 features.
// GLM result must equal the sum of per-instance multinomial_logit_lpmf calls.
// -----------------------------------------------------------------------
TEST(ProbMultinomialLogitGLM, K3_N3_manual) {
  Eigen::MatrixXd x(3, 2);
  x <<  1.0,  0.5,
       -0.5,  1.0,
        0.0, -1.0;

  Eigen::RowVectorXd alpha(3);
  alpha << 0.2, -0.1, 0.5;

  Eigen::MatrixXd beta(2, 3);
  beta <<  0.3, -0.2,  0.1,
          -0.1,  0.4, -0.3;

  std::vector<std::vector<int>> y{
      {2, 1, 3},
      {0, 4, 1},
      {3, 0, 2}};

  double expected = 0;
  for (int n = 0; n < 3; ++n) {
    Eigen::VectorXd eta = alpha.transpose() + beta.transpose() * x.row(n).transpose();
    expected += stan::math::multinomial_logit_lpmf(y[n], eta);
  }

  EXPECT_FLOAT_EQ(expected,
                  stan::math::multinomial_logit_glm_lpmf(y, x, alpha, beta));
}

// -----------------------------------------------------------------------
// Matrix alpha (N x K): per-instance intercepts, full design matrix.
// Result must equal the sum of per-instance multinomial_logit_lpmf calls
// each using its own alpha row.
// -----------------------------------------------------------------------
TEST(ProbMultinomialLogitGLM, matrixAlpha_fullX) {
  const int N = 4, K = 3, M = 2;
  Eigen::MatrixXd x(N, M);
  x <<  1.0,  0.5,
       -0.5,  1.0,
        0.0, -1.0,
        0.7,  0.3;

  Eigen::MatrixXd alpha(N, K);
  alpha <<  0.2, -0.1,  0.5,
           -0.3,  0.4,  0.1,
            0.1,  0.0, -0.2,
            0.5, -0.3,  0.2;

  Eigen::MatrixXd beta(M, K);
  beta <<  0.3, -0.2,  0.1,
          -0.1,  0.4, -0.3;

  std::vector<std::vector<int>> y{
      {2, 1, 3},
      {0, 4, 1},
      {3, 0, 2},
      {1, 2, 1}};

  double expected = 0;
  for (int n = 0; n < N; ++n) {
    Eigen::VectorXd eta = alpha.row(n).transpose()
                          + beta.transpose() * x.row(n).transpose();
    expected += stan::math::multinomial_logit_lpmf(y[n], eta);
  }

  EXPECT_FLOAT_EQ(expected,
                  stan::math::multinomial_logit_glm_lpmf(y, x, alpha, beta));
}

// -----------------------------------------------------------------------
// propto=true: with non-autodiff parameters the entire log PMF is constant
// w.r.t. parameters, so Stan math's early-return yields 0.
// propto=false: full log PMF must be finite and negative.
// -----------------------------------------------------------------------
TEST(ProbMultinomialLogitGLM, propto) {
  Eigen::MatrixXd x(2, 2);
  x << 1.0, 0.0, 0.0, 1.0;
  Eigen::RowVectorXd alpha(3);
  alpha << 0.1, 0.2, 0.3;
  Eigen::MatrixXd beta(2, 3);
  beta << 0.1, 0.2, 0.3, 0.4, 0.5, 0.6;

  std::vector<std::vector<int>> y{{1, 2, 3}, {2, 1, 0}};

  const double logp_full
      = stan::math::multinomial_logit_glm_lpmf<false>(y, x, alpha, beta);
  EXPECT_TRUE(std::isfinite(logp_full));
  EXPECT_LT(logp_full, 0.0);

  // With all-double parameters and propto=true, Stan math elides all terms
  // (the full log PMF is a constant) and returns 0.
  EXPECT_EQ(stan::math::multinomial_logit_glm_lpmf<true>(y, x, alpha, beta),
            0.0);
}

// -----------------------------------------------------------------------
// alpha[n,k] = -inf forces softmax probability to 0 for class k on instance n.
// When y[n,k] = 0 the 0*log(0) term must evaluate to 0 (not NaN).
// -----------------------------------------------------------------------
TEST(ProbMultinomialLogitGLM, negInfAlpha) {
  // N=2, K=3, M=2. Some alpha entries are -inf, forcing p=0 for those classes.
  const int N = 2, K = 3, M = 2;

  Eigen::MatrixXd x(N, M);
  x << 1.0,  0.5,
       0.3, -0.7;

  Eigen::MatrixXd beta(M, K);
  beta <<  0.3, -0.2,  0.1,
          -0.1,  0.4, -0.3;

  // Per-instance (N x K) alpha with -inf entries; y is 0 for those classes.
  Eigen::MatrixXd alpha(N, K);
  alpha <<  0.2, -0.1, -stan::math::INFTY,
            -stan::math::INFTY,  0.4,  0.1;

  std::vector<std::vector<int>> y{
      {2, 1, 0},
      {0, 3, 2}};

  const double logp = stan::math::multinomial_logit_glm_lpmf(y, x, alpha, beta);
  EXPECT_TRUE(std::isfinite(logp));
  EXPECT_LT(logp, 0.0);

  // Reference: -inf drops the class from the softmax, equivalent to
  // multinomial_logit_lpmf over the remaining classes.
  // lgamma(0 + 1) = 0 so the normalizing constant is unaffected.
  double expected = 0;
  for (int n = 0; n < N; ++n) {
    const Eigen::VectorXd eta
        = beta.transpose() * x.row(n).transpose() + alpha.row(n).transpose();
    std::vector<int> y_in;
    std::vector<double> eta_in;
    for (int k = 0; k < K; ++k) {
      if (std::isfinite(eta(k))) {
        y_in.push_back(y[n][k]);
        eta_in.push_back(eta(k));
      }
    }
    expected += stan::math::multinomial_logit_lpmf(
        y_in, Eigen::Map<const Eigen::VectorXd>(eta_in.data(), eta_in.size()));
  }

  EXPECT_FLOAT_EQ(expected, logp);
}

// -----------------------------------------------------------------------
// Error handling
// -----------------------------------------------------------------------
TEST(ProbMultinomialLogitGLM, throwsCorrectly) {
  Eigen::MatrixXd x(2, 2);
  x << 1.0, 0.5, -0.5, 1.0;
  Eigen::RowVectorXd alpha(3);
  alpha << 0.1, 0.2, 0.3;
  Eigen::MatrixXd beta(2, 3);
  beta << 0.1, 0.2, 0.3, 0.4, 0.5, 0.6;

  std::vector<std::vector<int>> y{{1, 2, 3}, {2, 1, 0}};

  // Mismatched number of instances (y has 1 row, x has 2)
  EXPECT_THROW(stan::math::multinomial_logit_glm_lpmf({{1, 2, 3}}, x, alpha,
                                                       beta),
               std::invalid_argument);

  // Mismatched number of classes in one outcome vector
  EXPECT_THROW(
      stan::math::multinomial_logit_glm_lpmf({{1, 2}, {2, 1, 0}}, x, alpha,
                                              beta),
      std::invalid_argument);

  // Negative counts
  EXPECT_THROW(
      stan::math::multinomial_logit_glm_lpmf({{-1, 2, 3}, {2, 1, 0}}, x,
                                              alpha, beta),
      std::domain_error);

  // Infinite alpha
  Eigen::RowVectorXd alpha_inf = alpha;
  alpha_inf[0] = stan::math::INFTY;
  EXPECT_THROW(
      stan::math::multinomial_logit_glm_lpmf(y, x, alpha_inf, beta),
      std::domain_error);

  // Infinite beta
  Eigen::MatrixXd beta_inf = beta;
  beta_inf(0, 0) = stan::math::INFTY;
  EXPECT_THROW(
      stan::math::multinomial_logit_glm_lpmf(y, x, alpha, beta_inf),
      std::domain_error);

  // Infinite x
  Eigen::MatrixXd x_inf = x;
  x_inf(0, 0) = stan::math::INFTY;
  EXPECT_THROW(
      stan::math::multinomial_logit_glm_lpmf(y, x_inf, alpha, beta),
      std::domain_error);

  // alpha cols mismatch with K (beta has 3 cols but alpha has 2 cols)
  Eigen::RowVectorXd alpha_bad(2);
  alpha_bad << 0.1, 0.2;
  EXPECT_THROW(
      stan::math::multinomial_logit_glm_lpmf(y, x, alpha_bad, beta),
      std::invalid_argument);

  // Matrix alpha with wrong number of rows
  Eigen::MatrixXd alpha_mat_bad_rows(3, 3);
  EXPECT_THROW(
      stan::math::multinomial_logit_glm_lpmf(y, x, alpha_mat_bad_rows, beta),
      std::invalid_argument);

  // Matrix alpha with wrong number of columns
  Eigen::MatrixXd alpha_mat_bad_cols(2, 2);
  EXPECT_THROW(
      stan::math::multinomial_logit_glm_lpmf(y, x, alpha_mat_bad_cols, beta),
      std::invalid_argument);

  // beta rows mismatch with x cols (x has 2 cols, beta has 3 rows)
  EXPECT_THROW(
      stan::math::multinomial_logit_glm_lpmf(y, x, alpha,
                                              Eigen::MatrixXd(3, 3)),
      std::invalid_argument);
}
