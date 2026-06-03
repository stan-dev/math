#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>

TEST_F(AgradRev, mathMixScalFun_multinomial_logit_glm_lpmf) {
  auto f = [](const auto& y) {
    return [=](const auto& x, const auto& alpha, const auto& beta) {
      return stan::math::multinomial_logit_glm_lpmf(y, x, alpha, beta);
    };
  };

  std::vector<std::vector<int>> y2{{2, 1}, {0, 3}};
  Eigen::MatrixXd x(2, 2);
  x << 0.5, -0.3, 0.1, 0.8;
  Eigen::RowVectorXd alpha_row(2);
  alpha_row << 0.2, -0.1;
  Eigen::MatrixXd alpha_mat(2, 2);
  alpha_mat << 0.2, -0.1, 0.3, 0.4;
  Eigen::MatrixXd beta(2, 2);
  beta << 0.3, -0.2, -0.1, 0.4;

  // broadcast row-vector alpha
  stan::test::expect_ad(f(y2), x, alpha_row, beta);
  // per-instance matrix alpha
  stan::test::expect_ad(f(y2), x, alpha_mat, beta);
}

// Boundary conditions of the multinomial-logit GLM distribution, checked for
// value + forward/reverse higher-order derivatives via expect_ad.
TEST_F(AgradRev, mathMixScalFun_multinomial_logit_glm_lpmf_boundaries) {
  auto f = [](const auto& y) {
    return [=](const auto& x, const auto& alpha, const auto& beta) {
      return stan::math::multinomial_logit_glm_lpmf(y, x, alpha, beta);
    };
  };

  // Shared K=3 design used by several cases below.
  Eigen::MatrixXd x(2, 2);
  x << 0.5, -0.3, 0.1, 0.8;
  Eigen::RowVectorXd alpha_row(3);
  alpha_row << 0.2, -0.1, 0.3;
  Eigen::MatrixXd alpha_mat(2, 3);
  alpha_mat << 0.2, -0.1, 0.3, -0.3, 0.4, 0.1;
  Eigen::MatrixXd beta(2, 3);
  beta << 0.3, -0.2, 0.1, -0.1, 0.4, -0.3;

  // Instance with zero total counts (S_n = 0): that row must contribute 0 to
  // the log-PMF and 0 to every gradient (delta_n = y_n - S_n p_n = 0).
  {
    std::vector<std::vector<int>> y{{0, 0, 0}, {2, 1, 3}};
    stan::test::expect_ad(f(y), x, alpha_row, beta);
    stan::test::expect_ad(f(y), x, alpha_mat, beta);
  }

  // Zero counts in individual classes: exercises lmultiply's 0*log(0)=0
  // handling under autodiff (the gradient contribution of a 0-count class is 0
  // regardless of its probability).
  {
    std::vector<std::vector<int>> y{{0, 5, 0}, {4, 0, 0}};
    stan::test::expect_ad(f(y), x, alpha_row, beta);
    stan::test::expect_ad(f(y), x, alpha_mat, beta);
  }

  // Single class (K = 1): degenerate simplex, softmax == 1 so the log-PMF and
  // all gradients are identically 0.
  {
    std::vector<std::vector<int>> y{{3}, {2}};
    Eigen::RowVectorXd alpha_k1(1);
    alpha_k1 << 0.2;
    Eigen::MatrixXd beta_k1(2, 1);
    beta_k1 << 0.3, -0.1;
    stan::test::expect_ad(f(y), x, alpha_k1, beta_k1);
  }

  // Single instance (N = 1).
  {
    std::vector<std::vector<int>> y{{2, 1, 3}};
    Eigen::MatrixXd x1(1, 2);
    x1 << 0.5, -0.3;
    stan::test::expect_ad(f(y), x1, alpha_row, beta);
  }

  // Large-magnitude linear predictor: drives some softmax probabilities near
  // 0/1, exercising the row-max numerical-stability shift (eta - rowwise max)
  // through the higher-order AD stack. Counts are kept off the vanishing-
  // probability classes so the reference log-PMF stays well conditioned.
  {
    std::vector<std::vector<int>> y{{3, 0, 4}, {0, 4, 1}};
    Eigen::MatrixXd x_big(2, 2);
    x_big << 2.5, -2.0, 1.8, 2.2;
    Eigen::RowVectorXd alpha_big_row(3);
    alpha_big_row << 1.5, -1.0, 2.0;
    Eigen::MatrixXd alpha_big_mat(2, 3);
    alpha_big_mat << 1.5, -1.0, 2.0, -1.2, 1.1, 0.7;
    Eigen::MatrixXd beta_big(2, 3);
    beta_big << 1.2, -1.5, 0.8, -1.0, 1.6, -1.3;
    stan::test::expect_ad(f(y), x_big, alpha_big_row, beta_big);
    stan::test::expect_ad(f(y), x_big, alpha_big_mat, beta_big);
  }
}
