#include <test/unit/math/test_ad.hpp>

TEST(mathMixScalFun, hypergeometric_1f0) {
  auto f = [](const auto& x1, const auto& x2) {
    using stan::math::hypergeometric_1f0;
    return hypergeometric_1f0(x1, x2);
  };

  stan::test::expect_ad(f, 5, 0.3);
  stan::test::expect_ad(f, 3.4, 0.9);
  stan::test::expect_ad(f, 3.4, 0.1);
  stan::test::expect_ad(f, 5, -0.7);
  stan::test::expect_ad(f, 7, -0.1);
  stan::test::expect_ad(f, 2.8, 0.8);
}
