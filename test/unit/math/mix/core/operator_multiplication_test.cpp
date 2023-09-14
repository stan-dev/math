#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/util.hpp>

TEST_F(mathMix, operatorMultiplication) {
  auto f = [](const auto& x1, const auto& x2) { return x1 * x2; };
  bool disable_lhs_int = true;
  stan::test::expect_common_binary(f, disable_lhs_int);
  stan::test::expect_complex_common_binary(f);
}
