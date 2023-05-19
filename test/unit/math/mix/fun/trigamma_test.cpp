#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, trigamma) {
  using stan::math::trigamma;
  using stan::test::ad_tolerances;
  using stan::test::expect_unary_vectorized;
  using stan::test::relative_tolerance;

  auto f = [](const auto& x1) { return trigamma(x1); };

  // reduce the tol_min for second and third order tests one order of magnitude
  stan::test::ad_tolerances tols;
  tols.hessian_hessian_ = relative_tolerance(1e-4, 1e-2);
  tols.hessian_fvar_hessian_ = relative_tolerance(1e-4, 1e-2);
  tols.grad_hessian_hessian_ = relative_tolerance(1e-3, 1e-2);
  tols.grad_hessian_grad_hessian_ = relative_tolerance(1e-2, 1e-1);

  expect_unary_vectorized<stan::test::ScalarSupport::Real>(
      tols, f, -103.52, -0.9, -0.5, 0, 0.5, 1.3, 5.1, 19.2);

  // reduce tol_min for first deriv tests one order, second derivs four orders,
  // and third derivs three orders
  stan::test::ad_tolerances tols2;
  tols2.gradient_grad_ = relative_tolerance(1e-4, 1e-3);
  tols2.gradient_fvar_grad_ = relative_tolerance(1e-4, 1e-3);
  tols2.hessian_grad_ = relative_tolerance(1e-4, 1e-3);
  tols2.hessian_fvar_grad_ = relative_tolerance(1e-4, 1e-3);
  tols2.hessian_hessian_ = relative_tolerance(1e2, 1);
  tols2.hessian_fvar_hessian_ = relative_tolerance(1e2, 1);
  tols2.grad_hessian_hessian_ = relative_tolerance(1, 1);
  tols2.grad_hessian_grad_hessian_ = relative_tolerance(1, 1);

  expect_unary_vectorized(tols2, f, -20, 1, 5);
}
