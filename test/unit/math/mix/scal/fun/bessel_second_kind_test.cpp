#include <test/unit/math/test_ad.hpp>
#include <limits>

TEST(mathMixScalFun, besselSecondKind) {
  // bind integer arg because can't autodiff through
  auto f = [](const int x1) {
    return
        [=](const auto& x2) { return stan::math::bessel_second_kind(x1, x2); };
  };
  stan::test::expect_ad(f(0), 4.0);
  stan::test::expect_ad(f(0), std::numeric_limits<double>::quiet_NaN());
  stan::test::expect_ad(f(1), -3.0);
  stan::test::expect_ad(f(1), 3.0);
  stan::test::expect_ad(f(1), std::numeric_limits<double>::quiet_NaN());
  stan::test::expect_ad(f(2), 2.79);
}
