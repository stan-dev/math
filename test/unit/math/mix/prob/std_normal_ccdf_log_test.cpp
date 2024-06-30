#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>

TEST_F(AgradRev, mathMixScalFun_std_normal_lccdf) {
  auto f = [](const auto& y) { return stan::math::std_normal_lccdf(y); };

  stan::test::expect_ad(f, -50.0);
  stan::test::expect_ad(f, -20.0 * stan::math::SQRT_TWO);
  stan::test::expect_ad(f, -5.5);
  stan::test::expect_ad(f, 0.0);
  stan::test::expect_ad(f, 0.15);
  stan::test::expect_ad(f, 1.14);
  stan::test::expect_ad(f, 3.00);
  stan::test::expect_ad(f, 10.00);
}
