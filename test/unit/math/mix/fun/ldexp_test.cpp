#include <test/unit/math/test_ad.hpp>
#include <stan/math/prim.hpp>

TEST(mathMixScalFun, ldexp) {
  auto f = [](const auto& x1) { return stan::math::ldexp(x1, 5); };

  stan::test::expect_ad(f, 3.1);
  stan::test::expect_ad(f, 0.0);
  stan::test::expect_ad(f, -1.5);
  stan::test::expect_ad(f, stan::math::INFTY);
  stan::test::expect_ad(f, stan::math::NOT_A_NUMBER);
}
