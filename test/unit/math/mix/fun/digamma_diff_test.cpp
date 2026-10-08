#include <test/unit/math/test_ad.hpp>

TEST(mathMixScalFun, digamma_diff) {
  auto f = [](const auto& x, const auto& d) {
    return stan::math::digamma_diff(x, d);
  };
  // both sides of the shift threshold x = 10, integer and real offsets
  stan::test::expect_ad(f, 0.5, 1.0);
  stan::test::expect_ad(f, 2.3, 0.5);
  stan::test::expect_ad(f, 9.5, 3.0);
  stan::test::expect_ad(f, 10.5, 0.25);
  stan::test::expect_ad(f, 25.0, 57.0);
  stan::test::expect_ad(f, 1e3, 117.0);
  // invalid arguments throw for every type
  stan::test::expect_ad(f, -1.0, 2.0);
  stan::test::expect_ad(f, 2.0, -1.0);
}
