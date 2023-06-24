#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>

TEST(mathMixScalFun, categorical_logit_glm_lpmf) {
  auto f = [](const int y) {
    return [=](const auto& x, const auto& alpha, const auto& beta) {
      return stan::math::categorical_logit_glm_lpmf(y, x, alpha, beta);
    };
  };

  Eigen::MatrixXd x = Eigen::MatrixXd::Random(2, 2);
  Eigen::VectorXd alpha = Eigen::VectorXd::Random(2);
  Eigen::VectorXd beta = Eigen::VectorXd::Random(2);

  stan::test::expect_ad(f(0), x, alpha, beta);
  stan::test::expect_ad(f(1), x, alpha, beta);
}
