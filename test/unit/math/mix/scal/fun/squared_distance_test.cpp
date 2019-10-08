#include <test/unit/math/test_ad.hpp>
#include <limits>

TEST(mathMixScalFun, squaredDistance) {
  auto f = [](const auto& x1, const auto& x2) {
    return stan::math::squared_distance(x1, x2);
  };
  stan::test::expect_ad(f, 1.0, 1.0);
  stan::test::expect_ad(f, 1.0, 4.0);
  stan::test::expect_ad(f, 4.0, 1.0);
  stan::test::expect_ad(f, 4.0, 4.0);

  double nan = std::numeric_limits<double>::quiet_NaN();
  stan::test::expect_ad(f, 1.0, nan);
  stan::test::expect_ad(f, nan, 1.0);
  stan::test::expect_ad(f, nan, nan);
}
