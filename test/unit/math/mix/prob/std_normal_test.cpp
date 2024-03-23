#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>

TEST_F(AgradRev, mathMixScalFun_std_normal) {
  auto f = [](const auto& y) { return stan::math::std_normal_lpdf(y); };

  stan::test::expect_ad(f, -0.3);
  stan::test::expect_ad(f, 0.0);
  stan::test::expect_ad(f, 1.7);
}
