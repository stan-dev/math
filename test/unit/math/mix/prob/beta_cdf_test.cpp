#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>

TEST_F(AgradRev, mathMixScalFun_beta_cdf) {
  auto f = [](const auto& y, const auto& alpha, const auto& beta) {
    return stan::math::beta_cdf(y, alpha, beta);
  };

  stan::test::expect_ad(f, 0.8, 1.1, 2.3);
  stan::test::expect_ad(f, 0.8, 1.2, 2.3);
  stan::test::expect_ad(f, 0.2, 6.1, 0.3);
  stan::test::expect_ad(f, 0.3, 12.1, 15.1);
}
