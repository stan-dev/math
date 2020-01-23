#include <test/unit/math/test_ad.hpp>
#include <limits>

TEST(mathMixScalFun, step) {
  auto f = [](const auto& x) { return stan::math::step(x); };
  stan::test::expect_ad(f, -18765.3);
  stan::test::expect_ad(f, -1.0);
  stan::test::expect_ad(f, -0.001);
  stan::test::expect_ad(f, 0.001);
  stan::test::expect_ad(f, 1.0);
  stan::test::expect_ad(f, 3.5);
  stan::test::expect_ad(f, 2713.5);

  stan::test::expect_ad(f, std::numeric_limits<double>::quiet_NaN());
}
