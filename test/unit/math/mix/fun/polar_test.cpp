#include <test/unit/math/test_ad.hpp>
#include <complex>
#include <iostream>
#include <vector>

TEST(mathMixFun, polar) {
  auto f = [](const auto& r, const auto& theta) {
    using stan::math::polar;
    auto y = polar(r, theta);
    return y;
  };

  stan::test::expect_ad(f, 0.2, 0.5);
  stan::test::expect_ad(f, 1, 0.5);
  stan::test::expect_ad(f, 0.2, 1);
  stan::test::expect_ad(f, 1, 1);

  auto f1 = [](const auto& r, const auto& theta) {
    using stan::math::polar;
    auto y = polar(r, theta);
    return y.real();
  };

  auto f2 = [](const auto& r, const auto& theta) {
    using stan::math::polar;
    auto y = polar(r, theta);
    return y.imag();
  };
  stan::test::expect_ad(f1, 1, 1);
  stan::test::expect_ad(f2, 1, 1);

  std::vector<double> rs{-2.3, -0.3, 0.4, 1.2};
  std::vector<double> thetas{-0.2, 0, 0.1};

  for (double r : rs) {
    for (double theta : thetas) {
      stan::test::expect_ad(f1, r, theta);
      stan::test::expect_ad(f2, r, theta);
    }
  }

  for (double r : rs) {
    stan::test::expect_ad(f1, r, 1);
    stan::test::expect_ad(f2, r, 1);
  }
  for (double theta : thetas) {
    stan::test::expect_ad(f1, 1, theta);
    stan::test::expect_ad(f2, 1, theta);
  }
}
