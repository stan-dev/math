#include <test/unit/math/test_ad.hpp>
#include <limits>

TEST(mathMixScalFun, binaryLogLoss) {
  // bind integer arg because can't autodiff
  auto f = [](const auto& x1) {
    return [=](const auto& x2) { return stan::math::binary_log_loss(x1, x2); };
  };
  for (int y = 0; y <= 1; ++y) {
    stan::test::expect_ad(f(y), std::numeric_limits<double>::quiet_NaN());
    for (double y_hat = 0.05; y_hat < 1.0; y_hat += 0.05) {
      std::vector<int> sti_y{y, y, y};
      std::vector<std::vector<int>> ststi_y{sti_y, sti_y, sti_y};
      Eigen::VectorXi ei_y = Eigen::VectorXi::Constant(3, y);
      std::vector<Eigen::VectorXi> stei_y{ei_y, ei_y, ei_y};
      Eigen::MatrixXd ed_yhat = Eigen::MatrixXd::Constant(3, 3, y_hat);
      std::vector<Eigen::VectorXd> sted_y{
          ed_yhat.diagonal(), ed_yhat.diagonal(), ed_yhat.diagonal()};
      stan::test::expect_ad(f(y), y_hat);
      stan::test::expect_ad(f(ei_y), y_hat);
      stan::test::expect_ad(f(ei_y), ed_yhat.diagonal().eval());
      stan::test::expect_ad(f(Eigen::MatrixXi::Constant(3, 3, y)), y_hat);
      stan::test::expect_ad(f(y), ed_yhat);
      stan::test::expect_ad(f(sti_y), y_hat);
      stan::test::expect_ad(f(ststi_y), y_hat);
      stan::test::expect_ad(f(stei_y), y_hat);
      stan::test::expect_ad(f(y), sted_y);
      stan::test::expect_ad(f(stei_y), sted_y);
    }
  }
}
