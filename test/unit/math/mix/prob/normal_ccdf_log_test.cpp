#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>

TEST_F(AgradRev, mathMixScalFun_normal_lccdf) {
  auto f = [](const double mu, const double sigma) {
    return
        [=](const auto& y) { return stan::math::normal_lccdf(y, mu, sigma); };
  };

  stan::test::expect_ad(f(0.0, 1.0), -50.0);
  stan::test::expect_ad(f(0.0, 1.0), -20.0 * stan::math::SQRT_TWO);
  stan::test::expect_ad(f(0.0, 1.0), -5.5);
  stan::test::expect_ad(f(0.0, 1.0), 0.0);
  stan::test::expect_ad(f(0.0, 1.0), 0.15);
  stan::test::expect_ad(f(0.0, 1.0), 1.14);
  stan::test::expect_ad(f(0.0, 1.0), 3.00);
  stan::test::expect_ad(f(0.0, 1.0), 10.00);
  stan::test::expect_ad(f(-1.0, 2.0), 1.50);
  stan::test::expect_ad(f(2.0, 1.0), 0.50);

  // third order autodiff tests can fail at borders of piecewise function
  stan::test::ad_tolerances tols;
  tols.grad_hessian_grad_hessian_ = 1e1;
  stan::test::expect_ad(tols, f(0.0, 1.0), 0.1 * stan::math::SQRT_TWO);
}
