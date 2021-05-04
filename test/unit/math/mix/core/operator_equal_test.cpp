#include <test/unit/math/test_ad.hpp>

TEST(mathMixCore, operatorEqual) {
  auto f = [](const auto& x1, const auto& x2) { return x1 == x2; };
  stan::test::expect_common_comparison(f);
  stan::test::expect_complex_common_comparison(f);
}
