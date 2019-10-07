#include <test/unit/math/test_ad.hpp>

TEST(mathMixScalFun, binaryLogLoss) {
  // bind integer arg because can't autodiff
  auto f = [](int x1) {
    return [&](const auto& x2) { return stan::math::binary_log_loss(x1, x2); };
  };
  for (int y = 0; y <= 1; ++y) {
    for (double y_hat = 0; y_hat <= 1.0; y_hat += 0.1) {
      stan::test::expect_ad(f(y), y_hat);
      stan::test::expect_ad(f(y), std::numeric_limits<double>::quiet_NaN());
    }
  }
}
