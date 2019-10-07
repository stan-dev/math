#include <test/unit/math/test_ad.hpp>

TEST(mathMixScalFun, fallingFactorial) {
  auto f = [](const int x2) {
    return
        [&](const auto& x1) { return stan::math::falling_factorial(x1, x2); };
  };
  stan::test::expect_ad(f(2), -3.0);  // throws
  stan::test::expect_ad(f(2), 4.0);
  stan::test::expect_ad(f(4), 4.0);
  stan::test::expect_ad(f(3), 5);
  stan::test::expect_ad(f(3), 5.0);
  stan::test::expect_ad(f(2), std::numeric_limits<double>::quiet_NaN());
}
