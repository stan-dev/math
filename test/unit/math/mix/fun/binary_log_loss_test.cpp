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
