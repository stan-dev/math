#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>

TEST(mathMixScalFun, std_normal_lcdf) {
  auto f = [](const auto& y) { return stan::math::std_normal_lcdf(y); };

  stan::test::expect_ad(f, -50.0);
  stan::test::expect_ad(f, -20.0 * stan::math::SQRT_TWO);
  stan::test::expect_ad(f, -5.5);
  stan::test::expect_ad(f, 0.0);
  stan::test::expect_ad(f, 0.15);
  stan::test::expect_ad(f, 1.14);
  stan::test::expect_ad(f, 3.00);
  stan::test::expect_ad(f, 10.00);

  // thid order autodiff tests can fail at borders of piecewise function
  stan::test::ad_tolerances tols;
  tols.grad_hessian_grad_hessian_ = 1e1;
  stan::test::expect_ad(tols, f, 0.1 * stan::math::SQRT_TWO);
}
