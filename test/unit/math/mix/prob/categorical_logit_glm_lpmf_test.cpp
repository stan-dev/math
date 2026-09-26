#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>

TEST_F(AgradRev, mathMixScalFun_categorical_logit_glm_lpmf) {
  auto f = [](const auto y) {
    return [=](const auto& x, const auto& alpha, const auto& beta) {
      return stan::math::categorical_logit_glm_lpmf(y, x, alpha, beta);
    };
  };

  std::vector<int> y{0, 1};
  Eigen::MatrixXd x = Eigen::MatrixXd::Random(2, 2);
  Eigen::RowVectorXd x_rowvec = x.row(0);
  Eigen::VectorXd alpha = Eigen::VectorXd::Random(2);
  Eigen::VectorXd beta = Eigen::VectorXd::Random(2);

  stan::test::expect_ad(f(y[0]), x, alpha, beta);
  stan::test::expect_ad(f(y[1]), x_rowvec, alpha, beta);
  stan::test::expect_ad(f(y), x, alpha, beta);
  stan::test::expect_ad(f(y), x_rowvec, alpha, beta);
}

// several categories: beta is a matrix, so the products inside are
// matrix-matrix products, which must also work with forward-mode scalars
TEST_F(AgradRev, mathMixScalFun_categorical_logit_glm_lpmf_matrix_beta) {
  auto f = [](const auto y) {
    return [=](const auto& x, const auto& alpha, const auto& beta) {
      return stan::math::categorical_logit_glm_lpmf(y, x, alpha, beta);
    };
  };

  std::vector<int> y{1, 3, 2};
  Eigen::MatrixXd x = Eigen::MatrixXd::Random(3, 2);
  Eigen::RowVectorXd x_rowvec = x.row(0);
  Eigen::VectorXd alpha = Eigen::VectorXd::Random(3);
  Eigen::MatrixXd beta = Eigen::MatrixXd::Random(2, 3);

  stan::test::expect_ad(f(y[0]), x_rowvec, alpha, beta);
  stan::test::expect_ad(f(y), x, alpha, beta);
  stan::test::expect_ad(f(y), x_rowvec, alpha, beta);
}
