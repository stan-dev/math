#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>

TEST(mathMixScalFun, normal_lpdf) {
  auto f = [](const double mu, const double sigma) {
    return [=](const auto& y) { return stan::math::normal_lpdf(y, mu, sigma); };
  };

  stan::test::expect_ad(f(0, 1), -2.3);
  stan::test::expect_ad(f(0, 1), 0.0);
  stan::test::expect_ad(f(0, 1), 1.7);

  Eigen::VectorXd y(1);
  y << 1.7;
  stan::test::expect_ad(f(0, 1), y);
  stan::test::expect_ad_matvar(f(0, 1), y);
}
