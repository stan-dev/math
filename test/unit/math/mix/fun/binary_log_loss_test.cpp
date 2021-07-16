#include <test/unit/math/test_ad.hpp>
#include <limits>

TEST(mathMixScalFun, binaryLogLoss) {
  // bind integer arg because can't autodiff
  auto f = [](int x1) {
    return [=](const auto& x2) { return stan::math::binary_log_loss(x1, x2); };
  };
  for (int y = 0; y <= 1; ++y) {
    stan::test::expect_ad(f(y), std::numeric_limits<double>::quiet_NaN());
    for (double y_hat = 0.05; y_hat < 1.0; y_hat += 0.05)
      stan::test::expect_ad(f(y), y_hat);
  }
}

TEST(mathMixScalFun, binaryLogLossvec) {
  auto f = [](const auto& x1, const auto& x2) {
    using stan::math::binary_log_loss;
    return binary_log_loss(x1, x2);
  };

  std::vector<int> std_in1{3, 1};
  Eigen::VectorXd in2(2);
  in2 << 0.5, 3.4;
  stan::test::expect_ad_vectorized_binary(f, std_in1, in2);

  std::vector<std::vector<int>> std_std_in1{std_in1, std_in1};
  Eigen::MatrixXd mat_in2 = in2.replicate(1, 2);
  stan::test::expect_ad_vectorized_binary(f, std_std_in1, mat_in2);
}

TEST(mathMixScalFun, binaryLogLossMatVar) {
  auto f = [](const auto& x1, const auto& x2) {
    using stan::math::binary_log_loss;
    return binary_log_loss(x1, x2);
  };

  std::vector<int> std_in1{3, 1};
  Eigen::VectorXd in2(2);
  in2 << 0.5, 3.4;
  stan::test::expect_ad_matvar(f, std_in1, in2);
  stan::test::expect_ad_matvar(f, std_in1[0], in2);
}

TEST(mathMixScalFun, binaryLogLossMatVarVec) {
  auto f = [](const auto& x1, const auto& x2) {
    using stan::math::binary_log_loss;
    return binary_log_loss(x1, x2);
  };

  std::vector<int> std_in1{3, 1};
  Eigen::MatrixXd in2(2, 2);
  in2 << 0.5, 3.4, 0.5, 3.5;
  stan::test::expect_ad_vectorized_matvar(f, std_in1, in2);
}
