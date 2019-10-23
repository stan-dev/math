#include <test/unit/math/test_ad.hpp>

TEST(mathMixScalFun, gammaP) {
  auto f = [](const auto& x1, const auto& x2) {
    return stan::math::gamma_p(x1, x2);
  };
  stan::test::expect_common_nonzero_binary(f);
  stan::test::expect_ad(f, 0.5001, 1.0001);
}
