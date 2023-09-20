#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(ProbBinomialLogitGLM, matchesNonGLM) {
  using stan::math::binomial_logit_glm_lpmf;
  using stan::math::binomial_logit_lpmf;

  std::vector<int> n{1, 2};
  std::vector<int> N{5, 4};
  Eigen::MatrixXd x = Eigen::MatrixXd::Random(2, 2);
  Eigen::RowVectorXd x_row = x.row(0);
  Eigen::VectorXd alpha = Eigen::VectorXd::Random(2);
  Eigen::VectorXd beta = Eigen::VectorXd::Random(2);

  Eigen::VectorXd theta = alpha + x * beta;

  EXPECT_FLOAT_EQ(binomial_logit_lpmf(n, N, theta),
                  binomial_logit_glm_lpmf(n, N, x, alpha, beta));
  EXPECT_FLOAT_EQ(binomial_logit_lpmf(n[0], N, theta),
                  binomial_logit_glm_lpmf(n[0], N, x, alpha, beta));
  EXPECT_FLOAT_EQ(binomial_logit_lpmf(n, N[0], theta),
                  binomial_logit_glm_lpmf(n, N[0], x, alpha, beta));
  EXPECT_FLOAT_EQ(binomial_logit_lpmf(n[0], N[0], theta),
                  binomial_logit_glm_lpmf(n[0], N[0], x, alpha, beta));

  theta = (alpha[0] + (x * beta).array()).matrix();

  EXPECT_FLOAT_EQ(binomial_logit_lpmf(n, N, theta),
                  binomial_logit_glm_lpmf(n, N, x, alpha[0], beta));
  EXPECT_FLOAT_EQ(binomial_logit_lpmf(n[0], N, theta),
                  binomial_logit_glm_lpmf(n[0], N, x, alpha[0], beta));
  EXPECT_FLOAT_EQ(binomial_logit_lpmf(n, N[0], theta),
                  binomial_logit_glm_lpmf(n, N[0], x, alpha[0], beta));
  EXPECT_FLOAT_EQ(binomial_logit_lpmf(n[0], N[0], theta),
                  binomial_logit_glm_lpmf(n[0], N[0], x, alpha[0], beta));

  theta = (alpha.array() + (x_row * beta)(0, 0)).matrix();

  EXPECT_FLOAT_EQ(binomial_logit_lpmf(n, N, theta),
                  binomial_logit_glm_lpmf(n, N, x_row, alpha, beta));
  EXPECT_FLOAT_EQ(binomial_logit_lpmf(n[0], N, theta),
                  binomial_logit_glm_lpmf(n[0], N, x_row, alpha, beta));
  EXPECT_FLOAT_EQ(binomial_logit_lpmf(n, N[0], theta),
                  binomial_logit_glm_lpmf(n, N[0], x_row, alpha, beta));
  EXPECT_FLOAT_EQ(binomial_logit_lpmf(n[0], N[0], theta),
                  binomial_logit_glm_lpmf(n[0], N[0], x_row, alpha, beta));
}

TEST(ProbBinomialLogitGLM, throwsCorrectly) {
  using stan::math::binomial_logit_glm_lpmf;
  using stan::math::INFTY;

  std::vector<int> n{1, 2};
  std::vector<int> N{5, 4};
  Eigen::MatrixXd x = Eigen::MatrixXd::Random(2, 2);
  Eigen::VectorXd alpha = Eigen::VectorXd::Random(2);
  Eigen::VectorXd beta = Eigen::VectorXd::Random(2);

  std::vector<int> N_mismatch_size{5, 4, 10};
  EXPECT_THROW(binomial_logit_glm_lpmf(n, N_mismatch_size, x, alpha, beta),
               std::invalid_argument);
  EXPECT_THROW(binomial_logit_glm_lpmf(500, 1, x, alpha, beta),
               std::domain_error);
  EXPECT_THROW(binomial_logit_glm_lpmf(-10, N, x, alpha, beta),
               std::domain_error);
  EXPECT_THROW(binomial_logit_glm_lpmf(n, -10, x, alpha, beta),
               std::domain_error);

  Eigen::VectorXd alpha_inf = alpha;
  alpha[0] = INFTY;
  Eigen::VectorXd beta_inf = beta;
  beta[0] = INFTY;
  Eigen::MatrixXd x_inf = x;
  x(0, 0) = INFTY;

  EXPECT_THROW(binomial_logit_glm_lpmf(n, N, x_inf, alpha, beta),
               std::domain_error);
  EXPECT_THROW(binomial_logit_glm_lpmf(n, N, x, alpha_inf, beta),
               std::domain_error);
  EXPECT_THROW(binomial_logit_glm_lpmf(n, N, x, alpha, beta_inf),
               std::domain_error);
}
