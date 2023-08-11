#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>

TEST(mathMixScalFun, poisson_log_glm_lpmf) {
  auto f = [](const auto y) {
    return [=](const auto& x, const auto& alpha, const auto& beta) {
      return stan::math::poisson_log_glm_lpmf(y, x, alpha, beta);
    };
  };

  std::vector<int> y{0, 1, 2};
  Eigen::MatrixXd x = Eigen::MatrixXd::Random(2, 2);
  Eigen::RowVectorXd x_rowvec = x.row(0);
  Eigen::VectorXd alpha = Eigen::VectorXd::Random(2);
  Eigen::VectorXd beta = Eigen::VectorXd::Random(2);

  stan::test::expect_ad(f(y[0]), x, alpha, beta);
  stan::test::expect_ad(f(y[0]), x, alpha[0], beta);

  stan::test::expect_ad(f(y[0]), x_rowvec, alpha, beta);
  stan::test::expect_ad(f(y[0]), x_rowvec, alpha[0], beta);

  stan::test::expect_ad(f(y), x, alpha, beta);
  stan::test::expect_ad(f(y), x, alpha[0], beta);

  stan::test::expect_ad(f(y), x_rowvec, alpha, beta);
  stan::test::expect_ad(f(y), x_rowvec, alpha[0], beta);
}
