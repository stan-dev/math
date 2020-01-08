#include <test/unit/math/test_ad.hpp>
#include <limits>

TEST(mathMixScalFun, fallingFactorial) {
  auto f = [](const int x2) {
    return
        [=](const auto& x1) { return stan::math::falling_factorial(x1, x2); };
  };
  stan::test::expect_ad(f(-2), -3.0);  // throws

  stan::test::expect_ad(f(3), 5);

  // 3rd order derivatives near zero
  stan::test::ad_tolerances tols;
  tols.grad_hessian_grad_hessian_ = 3.0;

  stan::test::expect_ad(tols, f(2), 4.0);
  stan::test::expect_ad(tols, f(4), 4.0);
  stan::test::expect_ad(tols, f(3), 5.0);

  stan::test::expect_ad(f(2), std::numeric_limits<double>::quiet_NaN());
}
