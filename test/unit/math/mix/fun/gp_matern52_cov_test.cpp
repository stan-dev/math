#include <test/unit/math/test_ad.hpp>
#include <vector>

TEST(mathMixMatFun, gpMatern52Cov) {
  auto f = [](const auto& x, const auto& sigma, const auto& l) {
    return stan::math::gp_matern52_cov(x, sigma, l);
  };

  double sigma = 0.2;
  double l = 5;
  std::vector<double> x1{-2};
  std::vector<double> x2{-2, 1};
  std::vector<double> x3{-2, -1, -0.5};
  stan::test::expect_ad(f, x1, sigma, l);
  stan::test::expect_ad(f, x2, sigma, l);
  stan::test::expect_ad(f, x3, sigma, l);

  std::vector<double> x0{};
  stan::test::expect_ad(f, x0, sigma, l);
}
