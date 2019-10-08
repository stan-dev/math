#include <test/unit/math/test_ad.hpp>
#include <limits>

TEST(mathMixScalFun, logInvLogitDiff) {
  auto f = [](const auto& x1, const auto& x2) {
    return stan::math::log_inv_logit_diff(x1, x2);
  };
  stan::test::expect_ad(f, 0.5, -1.0);
  stan::test::expect_ad(f, 0.5, 0.0);
  stan::test::expect_ad(f, 1.2, 0.6);
  stan::test::expect_ad(f, 3.4, 0.9);

  double nan = std::numeric_limits<double>::quiet_NaN();
  stan::test::expect_ad(f, 1.0, nan);
  stan::test::expect_ad(f, nan, 1.0);
  stan::test::expect_ad(f, nan, nan);
}
