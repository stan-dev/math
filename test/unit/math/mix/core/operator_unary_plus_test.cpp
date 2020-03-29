#include <test/unit/math/test_ad.hpp>

TEST(mathMixCore, operatorUnaryPlus) {
  auto f = [](const auto& x1) { return +x1; };
  stan::test::expect_common_unary(f);
  stan::test::expect_complex_common(f);
}
