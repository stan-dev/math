#include <test/unit/math/test_ad.hpp>

TEST(mathMixScalFun, logSumExp) {
  auto f = [](const auto& x1, const auto& x2) {
    return stan::math::log_sum_exp(x1, x2);
  };
  stan::test::expect_ad(f, 0.0, 0.0);
  stan::test::expect_ad(f, 0.5, -1.0);
  stan::test::expect_ad(f, 0.5, 0.0);
  stan::test::expect_ad(f, 0.5, 1.2);
  stan::test::expect_ad(f, 0.5, 1.4);
  stan::test::expect_ad(f, 1.0, 2.0);
  stan::test::expect_ad(f, 1.4, 0.5);
  stan::test::expect_ad(f, 1.2, 0.6);
  stan::test::expect_ad(f, 2.0, 1.0);
  stan::test::expect_ad(f, 3.0, 5.0);
  stan::test::expect_ad(f, 3.0, 6.0);
  stan::test::expect_ad(f, 3.4, 0.9);
  stan::test::expect_ad(f, 5.0, 2.0);

  stan::test::ad_tolerances tols;
  tols.gradient_grad_ = 0.01;
  tols.gradient_fvar_grad_ = 0.01;
  tols.hessian_grad_ = 0.01;
  tols.hessian_hessian_ = 3.0;
  tols.hessian_fvar_grad_ = 0.01;
  tols.hessian_fvar_hessian_ = 3.0;
  tols.grad_hessian_hessian_ = 3.0;
  tols.grad_hessian_grad_hessian_ = 3.0;
  stan::test::expect_ad(tols, f, -2.0, 15.0);
  stan::test::expect_ad(tols, f, 2.0, -15.0);

  double nan = std::numeric_limits<double>::quiet_NaN();
  stan::test::expect_ad(f, 1.0, nan);
  stan::test::expect_ad(f, nan, 1.0);
  stan::test::expect_ad(f, nan, nan);

  stan::test::expect_value(f, 1000.0, 10.0);
}
