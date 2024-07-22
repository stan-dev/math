#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>

TEST_F(AgradRev, mathMixScalFun_gamma_lcdf) {
  auto f = [](const auto& y, const auto& alpha, const auto& beta) {
    return stan::math::gamma_lcdf(y, alpha, beta);
  };

  stan::test::expect_ad(f, 0.8, 1.1, 2.3);
  stan::test::expect_ad(f, 0.8, 12, 2.3);
  stan::test::expect_ad(f, 5, 12, 2.3);
  stan::test::expect_ad(f, 5, 12, 15);
}
