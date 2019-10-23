#include <test/unit/math/test_ad.hpp>
#include <limits>

TEST(mathMixScalFun, lmgamma) {
  auto f = [](int x1) {
    return [=](const auto& x2) { return stan::math::hypot(x1, x2); };
  };
  stan::test::expect_ad(f(3), 3.2);
  stan::test::expect_ad(f(3), std::numeric_limits<double>::quiet_NaN());
}
