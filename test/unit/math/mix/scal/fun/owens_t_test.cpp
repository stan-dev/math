#include <test/unit/math/test_ad.hpp>
#include <limits>

TEST(mathMixScalFun, owensT) {
  auto f = [](const auto& x1, const auto& x2) {
    return stan::math::owens_t(x1, x2);
  };
  stan::test::expect_ad(f, -10.9, 13);
  stan::test::expect_ad(f, 0.5, 1.0);
  stan::test::expect_ad(f, 1.0, 2.0);

  double nan = std::numeric_limits<double>::quiet_NaN();
  stan::test::expect_ad(f, 1.0, nan);
  stan::test::expect_ad(f, nan, 1.0);
  stan::test::expect_ad(f, nan, nan);
}
