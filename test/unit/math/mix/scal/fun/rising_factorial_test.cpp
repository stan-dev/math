#include <test/unit/math/test_ad.hpp>
#include <limits>

TEST(mathMixScalFun, risingFactorial) {
  auto f = [](const int x2) {
    return [=](const auto& x1) { return stan::math::rising_factorial(x1, x2); };
  };

  stan::test::expect_ad(f(1), 4.0);
  stan::test::expect_ad(f(3), 5.0);
  stan::test::expect_ad(f(4), 4.0);

  // 3rd derivatives close to zero and rel tolerance fails
  stan::test::ad_tolerances tols;
  tols.grad_hessian_grad_hessian_ = 3.0;

  stan::test::expect_ad(tols, f(1), 5.0);
  stan::test::expect_ad(tols, f(2), 4.0);
  stan::test::expect_ad(tols, f(2), 4.0);
  stan::test::expect_ad(tols, f(4), 4.0);

  stan::test::expect_ad(f(2), std::numeric_limits<double>::quiet_NaN());
}
