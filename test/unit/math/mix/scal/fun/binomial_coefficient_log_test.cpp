#include <test/unit/math/test_ad.hpp>

TEST(mathMixScalFun, binomialCoefficientLog) {
  auto f = [](const auto& x1, const auto& x2) {
    return stan::math::binomial_coefficient_log(x1, x2);
  };
  stan::test::expect_ad(f, 3, 2);
  stan::test::expect_ad(f, 24.0, 12.0);
  stan::test::expect_common_nonzero_binary(f);
}
