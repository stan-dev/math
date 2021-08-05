#include <test/unit/math/test_ad.hpp>
#include <complex>

TEST(mathMixMatFun, to_complex) {
  auto f0 = []() { return stan::math::to_complex(); };
  auto f1 = [](const auto& x) { return stan::math::to_complex(x); };
  auto f2 = [](const auto& x, const auto& y) { return stan::math::to_complex(x, y); };

  //stan::test::expect_common_nullary(f0);
  //stan::test::expect_ad(f0);
  stan::test::expect_common_unary(f1);
  stan::test::expect_common_binary(f2);
  
}
