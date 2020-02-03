#include <test/unit/math/test_ad.hpp>
#include <limits>
#include <vector>

TEST(mathMixScalFun, logMix) {
  auto f = [](const auto& x1, const auto& x2, const auto& x3) {
    return stan::math::log_mix(x1, x2, x3);
  };
  stan::test::ad_tolerances tols;
  tols.hessian_hessian_ = 1e-2;
  tols.hessian_fvar_hessian_ = 1e-2;

  stan::test::expect_ad(tols, f, 0.0001, 0.1, 5.0);
  stan::test::expect_ad(tols, f, 0.1, -1.4, 3.99);
  stan::test::expect_ad(tols, f, 0.2, -2.0, 6.0);
  stan::test::expect_ad(tols, f, 0.2, 2.0, -6.0);
  stan::test::expect_ad(tols, f, 0.3, -10.4, 7.8);
  stan::test::expect_ad(tols, f, 0.3, 1.0, -2.0);
  stan::test::expect_ad(tols, f, 0.3, 1.0, 2.0);
  stan::test::expect_ad(tols, f, 0.3, 2.0, -6.0);
  stan::test::expect_ad(tols, f, 0.3, 2.0, 6.0);
  stan::test::expect_ad(tols, f, 0.3, 10.4, -10.9);
  stan::test::expect_ad(tols, f, 0.3, -1.4, 1.7);
  stan::test::expect_ad(tols, f, 0.6, 0.3, 0.5);
  stan::test::expect_ad(tols, f, 0.7, -4.7, 10.1);
  stan::test::expect_ad(tols, f, 0.7, -1.5, 2.0);
  stan::test::expect_ad(tols, f, 0.7, 0.1, 5.0);
  stan::test::expect_ad(tols, f, 0.7, 1.4, -1.9);
  stan::test::expect_ad(tols, f, 0.7, 1.4, 3.99);
  stan::test::expect_ad(tols, f, 0.7, 1.5, -2.0);
  stan::test::expect_ad(tols, f, 0.7, 1.5, 5.0);
  stan::test::expect_ad(tols, f, 0.7, 2.0, -6.0);
  stan::test::expect_ad(tols, f, 0.7, 2.0, 6.0);
  stan::test::expect_ad(tols, f, 0.7, 3.0, 5.0);
  stan::test::expect_ad(tols, f, 0.7, 4.0, 5.0);
  stan::test::expect_ad(tols, f, 0.999, 0.1, 5.0);

  // integer instantiations
  stan::test::expect_ad(tols, f, 0.5, 1, 1.0);
  stan::test::expect_ad(tols, f, 0.5, 1.0, 1);
  stan::test::expect_ad(tols, f, 0.5, 1, 1);

  // won't go out of bounds on finite diff past 0 or 1
  auto g = [](const auto& x1, const auto& x2, const auto& x3) {
    return stan::math::log_mix(x1 < 0 ? 0 : x1 >= 1 ? 1 : x1, x2, x3);
  };

  // requires special hessian and grad hessian tolerances;  grads ok
  stan::test::ad_tolerances tols2;
  tols2.hessian_hessian_ = 1e1;
  tols2.hessian_fvar_hessian_ = 1e1;
  tols2.grad_hessian_grad_hessian_ = 1e6;
  stan::test::expect_ad(tols2, g, 0, 1.0, 1.0);
  stan::test::expect_ad(tols2, g, 0, 1.0, 1);
  stan::test::expect_ad(tols2, g, 0, 1, 1.0);
  stan::test::expect_ad(tols2, g, 0, 1, 1);
}
