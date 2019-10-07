#include <test/unit/math/test_ad.hpp>

TEST(mixScalFun, abs) {
  auto f = [](const auto& x) { return stan::math::abs(x); };
  stan::test::expect_common_nonzero_unary(f);
  stan::test::expect_value(f, 0);
  stan::test::expect_value(f, 0.0);

  stan::test::expect_ad(f, -3);
  stan::test::expect_ad(f, -2);
  stan::test::expect_ad(f, 2);

  stan::test::expect_ad(f, -17.3);
  stan::test::expect_ad(f, -0.68);
  stan::test::expect_ad(f, 0.68);
  stan::test::expect_ad(f, 2.0);
  stan::test::expect_ad(f, 4.0);
}
