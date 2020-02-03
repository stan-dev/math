#include <test/unit/math/test_ad.hpp>
#include <limits>

TEST(mathMixScalFun, modifiedBesselSecondKind) {
  auto f = [](const int x1) {
    return [=](const auto& x2) {
      return stan::math::modified_bessel_second_kind(x1, x2);
    };
  };
  stan::test::expect_ad(f(1), -3.0);  // error

  stan::test::expect_ad(f(1), 4.0);
  stan::test::expect_ad(f(2), 2.0);

  stan::test::expect_ad(f(1), std::numeric_limits<double>::quiet_NaN());
}
