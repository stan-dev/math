#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>

TEST_F(AgradRev, mathMixScalFun_normal_id_glm_lpdf) {
  auto f = [](const auto& y, const auto& x) {
    return [=](const auto& alpha, const auto& beta, const auto& sigma) {
      return stan::math::normal_id_glm_lpdf(y, x, alpha, beta, sigma);
    };
  };
  auto f2 = [](const auto& beta, const auto& sigma) {
    return [=](const auto& y, const auto& x, const auto& alpha) {
      return stan::math::normal_id_glm_lpdf(y, x, alpha, beta, sigma);
    };
  };

  Eigen::VectorXd y = Eigen::VectorXd::Random(2);
  Eigen::MatrixXd x = Eigen::MatrixXd::Random(2, 2);
  Eigen::RowVectorXd x_rowvec = x.row(0);
  Eigen::VectorXd alpha = Eigen::VectorXd::Random(2);
  Eigen::VectorXd beta = Eigen::VectorXd::Random(2);
  Eigen::VectorXd sigma = Eigen::VectorXd::Random(2);

  stan::test::expect_ad(f(y, x), alpha, beta, sigma);
  stan::test::expect_ad(f(y, x), alpha[0], beta, sigma);
  stan::test::expect_ad(f(y, x), alpha, beta, sigma[0]);

  stan::test::expect_ad(f(y, x_rowvec), alpha, beta, sigma);
  stan::test::expect_ad(f(y, x_rowvec), alpha[0], beta, sigma);
  stan::test::expect_ad(f(y, x_rowvec), alpha, beta, sigma[0]);

  stan::test::expect_ad(f2(beta, sigma), y, x, alpha);
  stan::test::expect_ad(f2(beta, sigma[0]), y, x, alpha);
  stan::test::expect_ad(f2(beta, sigma), y, x, alpha[0]);

  stan::test::expect_ad(f2(beta, sigma), y, x_rowvec, alpha);
  stan::test::expect_ad(f2(beta, sigma[0]), y, x_rowvec, alpha);
  stan::test::expect_ad(f2(beta, sigma), y, x_rowvec, alpha[0]);
}
