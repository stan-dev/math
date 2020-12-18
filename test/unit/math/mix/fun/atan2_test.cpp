#include <test/unit/math/test_ad.hpp>
#include <limits>

TEST(mathMixCore, atan2) {
  auto f = [](const auto& x1, const auto& x2) {
    using std::atan2;
    using stan::math::atan2;
    return atan2(x1, x2);
  };
  double nan = std::numeric_limits<double>::quiet_NaN();
  stan::test::expect_ad(f, nan, 1.0);
  stan::test::expect_ad(f, 1.0, nan);
  stan::test::expect_ad(f, nan, nan);
  stan::test::expect_ad(f, 1.0, 1.0);
  stan::test::expect_ad(f, 1.5, 1.5);
  stan::test::expect_ad(f, 1.2, 3.9);
  stan::test::expect_ad(f, 0.5, 2.3);
}
