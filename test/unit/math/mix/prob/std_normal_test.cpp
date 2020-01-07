#include <test/unit/math/test_ad.hpp>
#include <limits>

TEST(mathMixScalFun, std_normal_derivatives) {
  auto f = [](const auto& y) { return stan::math::std_normal_lpdf(y); };

  stan::test::expect_ad(f, -0.3);
  stan::test::expect_ad(f, 0.0);
  stan::test::expect_ad(f, 1.7);
}
