#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>

TEST(mathMixScalFun, weibull_cdf) {
  auto f = [](const double alpha, const double sigma) {
    return
        [=](const auto& y) { return stan::math::weibull_cdf(y, alpha, sigma); };
  };

  stan::test::expect_ad(f(0.1, 1.0), 0.1);
  stan::test::expect_ad(f(1.0, 1.0), 0.001);
  stan::test::expect_ad(f(0.1, 1.0), 1.25);
  stan::test::expect_ad(f(1.0, 1.0), 0.15);
  stan::test::expect_ad(f(0.5, 0.5), 0.25);
  stan::test::expect_ad(f(1.0, 0.5), 2.0);
  stan::test::expect_ad(f(2.5, 3.0), 1.5);
}
