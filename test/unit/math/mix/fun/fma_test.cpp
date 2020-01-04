#include <test/unit/math/test_ad.hpp>
#include <limits>

TEST(mathMixScalFun, fma) {
  auto f = [](const auto& x1, const auto& x2, const auto& x3) {
    return stan::math::fma(x1, x2, x3);
  };

  // test all int instantiation patterns
  stan::test::expect_ad(f, 1.0, 2.0, 3.0);
  stan::test::expect_ad(f, 1.0, 2.0, 3);
  stan::test::expect_ad(f, 1.0, 2, 3.0);
  stan::test::expect_ad(f, 1.0, 2, 3);
  stan::test::expect_ad(f, 1, 2.0, 3.0);
  stan::test::expect_ad(f, 1, 2.0, 3);
  stan::test::expect_ad(f, 1, 2, 3.0);
  stan::test::expect_ad(f, 1, 2, 3);

  // replicate pre-existing test cases
  stan::test::expect_ad(f, 3.0, 5.0, 7.0);
  stan::test::expect_ad(f, 0.5, 1.2, 1.8);
  stan::test::expect_ad(f, 2.5, 1.5, 1.7);

  // test all nan instantiations
  double nan = std::numeric_limits<double>::quiet_NaN();
  stan::test::expect_ad(f, 2.5, 1.7, 1.5);
  stan::test::expect_ad(f, 2.5, 1.7, nan);
  stan::test::expect_ad(f, 2.5, nan, 1.5);
  stan::test::expect_ad(f, 2.5, nan, nan);
  stan::test::expect_ad(f, nan, 1.7, 1.5);
  stan::test::expect_ad(f, nan, 1.7, nan);
  stan::test::expect_ad(f, nan, nan, 1.5);
  stan::test::expect_ad(f, nan, nan, nan);
}
