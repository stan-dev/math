#include <test/unit/math/test_ad.hpp>
#include <limits>

TEST(mathMixScalFun, fmin) {
  auto f
      = [](const auto& x1, const auto& x2) { return stan::math::fmin(x1, x2); };
  stan::test::expect_ad(f, -3.0, 4.0);
  stan::test::expect_ad(f, 1.3, 2.0);
  stan::test::expect_ad(f, 2.0, -3.0);
  stan::test::expect_ad(f, 2.3, 2.0);
  stan::test::expect_ad(f, 2.5, 1.5);
  stan::test::expect_ad(f, 4.0, -3.0);

  double nan = std::numeric_limits<double>::quiet_NaN();
  stan::test::expect_ad(f, 1.0, nan);
  stan::test::expect_ad(f, nan, 1.0);
  stan::test::expect_ad(f, nan, nan);

  stan::test::expect_value(f, 2.0, 2.0);
}
