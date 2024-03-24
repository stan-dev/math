#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>

TEST(mathMixScalFun, weibull_cdf) {
  auto f = [](const auto& y, const auto& alpha, const auto& sigma) {
    return stan::math::weibull_cdf(y, alpha, sigma);
  };

  stan::test::expect_ad(f, 10, 10, 10);
}
