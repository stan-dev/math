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
