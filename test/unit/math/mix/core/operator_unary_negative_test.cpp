#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/util.hpp>

TEST_F(mathMix, operatorUnaryNegative) {
  auto f = [](const auto& x1) { return -x1; };
  stan::test::expect_common_unary(f);
  stan::test::expect_complex_common(f);
}
