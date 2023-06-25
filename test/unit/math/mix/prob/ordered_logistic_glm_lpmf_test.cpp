#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>

TEST(mathMixScalFun, ordered_logistic_glm_lpmf) {
  auto f = [](const int y) {
    return [=](const auto& x, const auto& beta, const auto& cutpoints) {
      return stan::math::ordered_logistic_glm_lpmf(y, x, beta, cutpoints);
    };
  };

  Eigen::MatrixXd x = Eigen::MatrixXd::Random(2, 2);
  Eigen::VectorXd cutpoints = Eigen::VectorXd::Random(2);
  Eigen::VectorXd beta = Eigen::VectorXd::Random(2);

  stan::test::expect_ad(f(1), x, beta, cutpoints);
  stan::test::expect_ad(f(2), x, beta, cutpoints);
}
